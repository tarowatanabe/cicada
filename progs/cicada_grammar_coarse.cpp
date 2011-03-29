//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// learn coarse grammar
// we assume that non-terminals are annotated by @id, and
// simply erasing bits will uncover reduced grammar
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>

#include <stdexcept>
#include <vector>
#include <deque>

#include <cicada/hypergraph.hpp>
#include <cicada/rule.hpp>

#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/math/special_functions/expm1.hpp>

#include <utils/bithack.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/compress_stream.hpp>
#include <utils/resource.hpp>
#include <utils/mathop.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/array_power2.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

typedef boost::filesystem::path path_type;

typedef cicada::HyperGraph hypergraph_type;
typedef hypergraph_type::rule_type     rule_type;
typedef hypergraph_type::rule_ptr_type rule_ptr_type;

typedef rule_type::symbol_type     symbol_type;
typedef rule_type::symbol_set_type symbol_set_type;

template <typename Tp>
struct ptr_hash : public boost::hash<Tp>
{
  typedef boost::hash<Tp> hasher_type;

  size_t operator()(const Tp* x) const
  {
    return (x ? hasher_type::operator()(*x) : size_t(0));
  }
  
  size_t operator()(const boost::shared_ptr<Tp>& x) const
  {
    return (x ? hasher_type::operator()(*x) : size_t(0));
  }

};

template <typename Tp>
struct ptr_equal
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x == y || (x && y && *x == *y);
  }
  
  bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
  {
    return x == y || (x && y && *x == *y);
  }
};

typedef double count_type;
typedef double prob_type;
typedef double logprob_type;

class Grammar : public google::dense_hash_map<rule_ptr_type, count_type, ptr_hash<rule_type>, ptr_equal<rule_type> >
{
public:
  typedef google::dense_hash_map<rule_ptr_type, count_type, ptr_hash<rule_type>, ptr_equal<rule_type> > count_set_type;
  
public:
  Grammar() : count_set_type() { count_set_type::set_empty_key(rule_ptr_type()); }
};

class Lexicon : public google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> >
{
public:
  typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> >lexicon_type;

  Lexicon() : lexicon_type() { lexicon_type::set_empty_key(symbol_type()); }
};

typedef Grammar grammar_type;
typedef Lexicon lexicon_type;

path_type input_grammar_file = "-";
path_type input_lexicon_file = "-";
path_type output_prefix;

symbol_type goal = "[ROOT]";

int max_order = 6;
int max_iteration = 25;

// naive variational bayes for smoothing... otherwise, dirichlet prior
bool variational_bayes_mode = false;
bool maximum_mode = false;

double prior = 0.01;

int threads = 1;

int debug = 0;

void options(int argc, char** argv);
void write_grammar(const path_type& path, const grammar_type& grammar);
void read_grammar(const path_type& path, grammar_type& grammar);
void read_lexicon(const path_type& path, lexicon_type& lexicon);

struct SimpleSymbol
{
  SimpleSymbol(const lexicon_type& __lexicon,
	       const symbol_type& __goal)
    : lexicon(__lexicon),
      goal(__goal) {}

  symbol_type operator()(const symbol_type& symbol)
  {
    if (! symbol.is_non_terminal()) return symbol;
    if (symbol == goal) return symbol;
    if (lexicon.find(symbol) != lexicon.end()) symbol;
    
    const utils::piece piece = symbol.non_terminal_strip();
    
    // default X
    return (piece.find('^') != utils::piece::npos() ? "[x^]" : "[x]");
  }
  
  const lexicon_type& lexicon;
  const symbol_type goal;
};

struct CoarseSymbol
{
  CoarseSymbol(const lexicon_type& __lexicon,
	       const symbol_type& __goal,
	       const int __bits)
    : lexicon(__lexicon),
      goal(__goal),
      bits(__bits) {}

  symbol_type operator()(const symbol_type& symbol)
  {
    if (! symbol.is_non_terminal()) return symbol;
    if (symbol == goal) return symbol;
    if (lexicon.find(symbol) != lexicon.end()) symbol;
    
    const size_t cache_pos = hash_value(symbol) & (caches.size() - 1);
    cache_type& cache = caches[cache_pos];
    if (cache.symbol != symbol) {
      namespace xpressive = boost::xpressive;
      
      typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
      typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
      
      static pregex re = (xpressive::s1= +(~xpressive::_s)) >> '@' >> (xpressive::s2= -+xpressive::_d);
      
      const utils::piece piece = symbol.non_terminal_strip();
      const int mask = 1 << bits;
      
      pmatch what;
      if (xpressive::regex_match(piece, what, re)) {
	const int value = (utils::lexical_cast<int>(what[2]) & (~mask));
	cache.annotated = '[' + what[1] + '@' + utils::lexical_cast<std::string>(value) + ']';
      } else
	cache.annotated = '[' + piece + "@0]";
      
      cache.symbol = symbol;
    }
    return cache.annotated;
  }
  
  struct Cache
  {
    symbol_type symbol;
    symbol_type annotated;
    
    Cache() : symbol(), annotated() {}
  };
  typedef Cache cache_type;
  typedef utils::array_power2<cache_type, 1024 * 8, std::allocator<cache_type> > cache_set_type;
  
  const lexicon_type& lexicon;
  const symbol_type goal;
  const int bits;
  
  cache_set_type caches;
};

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (output_prefix.empty())
      throw std::runtime_error("empty output prefix");
    
    grammar_type grammar;
    lexicon_type lexicon;
    
    read_grammar(input_grammar_file, grammar);
    
    read_lexicon(input_lexicon_file, lexicon);
    
    for (int order = max_order - 1; order >= 0; -- order) {
      
      // dump grammar...!
    }
    
    // finally reduce to minus-grammar:
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Tp>
struct greater_ptr_second
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x->second > y->second;
  }
};

template <typename Tp>
struct less_ptr_second
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x->second < y->second;
  }
};

void write_grammar(const path_type& path, const grammar_type& grammar)
{
#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<symbol_type, grammar_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
    std::allocator<std::pair<const symbol_type, grammar_type> > > count_set_type;
#else
  typedef sgi::hash_map<symbol_type, grammar_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
    std::allocator<std::pair<const symbol_type, grammar_type> > > count_set_type;
#endif
  
  typedef std::vector<const grammar_type::value_type*, std::allocator<const grammar_type::value_type*> > sorted_type;
  
  if (grammar.empty()) return;

  count_set_type counts;
  sorted_type sorted;
    
  grammar_type::const_iterator giter_end = grammar.end();
  for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
    counts[giter->first->lhs][giter->first] = giter->second;

  utils::compress_ostream os(path, 1024 * 1024);
  
  count_set_type::const_iterator citer_end = counts.end();
  for (count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
    const grammar_type& grammar = citer->second;
      
    sorted.clear();
      
    grammar_type::const_iterator giter_end = grammar.end();
    for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
      sorted.push_back(&(*giter));
      
    std::sort(sorted.begin(), sorted.end(), greater_ptr_second<grammar_type::value_type>());

    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      os << *((*siter)->first) << " ||| ||| " << (*siter)->second << '\n';
  }
}

inline
bool parse_rule(const std::string& line, rule_type& rule, double& logprob)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  namespace phoenix = boost::phoenix;

  std::string::const_iterator iter = line.begin();
  std::string::const_iterator iter_end = line.end();
  
  // parse rule-part, parse "|||" and "|||", then parse double
  
  return (rule.assign(iter, iter_end)
	  && qi::phrase_parse(iter, iter_end, qi::lit("|||") >> qi::lit("|||") >> qi::double_, standard::space, logprob)
	  && iter == iter_end);
}

void read_grammar(const path_type& path, grammar_type& grammar)
{
  utils::compress_istream is(path, 1024 * 1024);
  
  std::string line;
  rule_type rule;
  double logprob;
  
  while (std::getline(is, line)) {
    if (! parse_rule(line, rule, logprob)) continue;
    
    grammar[rule_type::create(rule)] = logprob;
  }
  
}

void read_lexicon(const path_type& path, lexicon_type& lexicon)
{
  utils::compress_istream is(path, 1024 * 1024);
  
  std::string line;
  rule_type rule;
  double logprob;
  
  while (std::getline(is, line)) {
    if (! parse_rule(line, rule, logprob)) continue;
    
    lexicon.insert(rule.lhs);
  }
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  
  desc.add_options()
    ("grammar", po::value<path_type>(&input_grammar_file)->default_value("-"), "input grammar")
    ("lexicon", po::value<path_type>(&input_lexicon_file)->default_value("-"), "input lexical rules")
    ("prefix",  po::value<path_type>(&output_prefix),      "output prefix")
    
    ("goal", po::value<symbol_type>(&goal)->default_value(goal), "goal")
    
    ("max-order",     po::value<int>(&max_order)->default_value(max_order),         "maximum order")
    ("max-iteration", po::value<int>(&max_iteration)->default_value(max_iteration), "maximum iterations")
    
    ("variational-bayes", po::bool_switch(&variational_bayes_mode), "variational Bayes estimates")
    ("maximum",           po::bool_switch(&maximum_mode),           "maximum estimates")
    
    ("prior",      po::value<double>(&prior)->default_value(prior), "Dirichlet prior")
    
    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::variables_map variables;
  
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}

