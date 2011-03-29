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
#include <cicada/semiring.hpp>

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

typedef cicada::semiring::Logprob<double> weight_type;

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

class ExpectedCounts : public google::dense_hash_map<symbol_type, weight_type, boost::hash<symbol_type>, std::equal_to<symbol_type> >
{
public:
  typedef google::dense_hash_map<symbol_type, weight_type, boost::hash<symbol_type>, std::equal_to<symbol_type> >expected_counts_type;
  
  ExpectedCounts() : expected_counts_type() { expected_counts_type::set_empty_key(symbol_type()); }
};


typedef Grammar grammar_type;
typedef Lexicon lexicon_type;
typedef ExpectedCounts expected_counts_type;

#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<symbol_type, grammar_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				  std::allocator<std::pair<const symbol_type, grammar_type> > > count_set_type;
#else
  typedef sgi::hash_map<symbol_type, grammar_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
			std::allocator<std::pair<const symbol_type, grammar_type> > > count_set_type;
#endif

path_type input_grammar_file = "-";
path_type input_lexicon_file = "-";
path_type output_file;

symbol_type goal = "[ROOT]";

int max_order = 6;
int max_iteration = 32;

// naive variational bayes for smoothing... otherwise, dirichlet prior
bool variational_bayes_mode = false;
bool maximum_mode = false;

double prior = 0.01;

int threads = 1;

int debug = 0;

void grammar_counts(const grammar_type& grammar, const lexicon_type& lexicon, expected_counts_type& counts);
template <typename Coarser>
void grammar_coarse(const grammar_type& grammar, const expected_counts_type& counts, grammar_type& coarse, Coarser coarser);

void write_grammar(const path_type& prefix, const int order, const grammar_type& grammar);
void read_grammar(const path_type& path, grammar_type& grammar);
void read_lexicon(const path_type& path, lexicon_type& lexicon);

void options(int argc, char** argv);

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
    if (lexicon.find(symbol) != lexicon.end()) return symbol;
    
    const utils::piece piece = symbol.non_terminal_strip();
    
    // default X
    return (piece.find('^') != utils::piece::npos() ? "[x^]" : "[x]");
  }
  
  const lexicon_type& lexicon;
  symbol_type goal;
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
    if (lexicon.find(symbol) != lexicon.end()) return symbol;
    
    const size_t cache_pos = hash_value(symbol) & (caches.size() - 1);
    cache_type& cache = caches[cache_pos];
    if (cache.symbol != symbol) {
      namespace xpressive = boost::xpressive;
      
      typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
      typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
      
      static pregex re = (xpressive::s1= +(~xpressive::_s)) >> '@' >> (xpressive::s2= -+xpressive::_d);
      
      const utils::piece piece = symbol.non_terminal_strip();
      const int mask = (1 << bits) - 1;
      
      pmatch what;
      if (xpressive::regex_match(piece, what, re)) {
	const int value = (utils::lexical_cast<int>(what[2]) & mask);
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
  symbol_type goal;
  int bits;
  
  cache_set_type caches;
};

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (output_file.empty())
      throw std::runtime_error("empty output file");

    threads = utils::bithack::max(threads, 1);
    
    grammar_type grammar;
    lexicon_type lexicon;
    
    read_grammar(input_grammar_file, grammar);
    
    read_lexicon(input_lexicon_file, lexicon);
    
    // compute expected counts over the grammar

    if (debug)
      std::cerr << "computing expected counts" << std::endl;
    
    expected_counts_type counts;
    grammar_counts(grammar, lexicon, counts);
    
    grammar_type coarse;
    for (int order = max_order - 1; order >= 0; -- order) {
      if (debug)
	std::cerr << "coarse order: " << order << std::endl;

      grammar_coarse(grammar, counts, coarse, CoarseSymbol(lexicon, goal, order));
      
      write_grammar(output_file, order + 1, coarse);
    }
    
    // finally reduce to minus-grammar:-)
    
    if (debug)
      std::cerr << "final coarse grammar" << std::endl;

    grammar_coarse(grammar, counts, coarse, SimpleSymbol(lexicon, goal));
    
    write_grammar(output_file, 0, coarse);
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Coarser>
void grammar_coarse(const grammar_type& grammar, const expected_counts_type& expected_counts, grammar_type& coarse, Coarser coarser)
{
  count_set_type counts;
  
  grammar_type::const_iterator riter_end = grammar.end();
  for (grammar_type::const_iterator riter = grammar.begin(); riter != riter_end; ++ riter) {
    
    // if not found... not reachable from ROOT!
    expected_counts_type::const_iterator eiter = expected_counts.find(riter->first->lhs);
    if (eiter == expected_counts.end())
      continue;
    
    const weight_type count = cicada::semiring::traits<weight_type>::exp(riter->second) * eiter->second;
    
    //
    // transform rule into a coarse rule
    //
    
    const symbol_type lhs = coarser(riter->first->lhs);
    symbol_set_type rhs(riter->first->rhs);
    symbol_set_type::iterator siter_end = rhs.end();
    for (symbol_set_type::iterator siter = rhs.begin(); siter != siter_end; ++ siter)
      *siter = coarser(*siter);
    
    std::pair<grammar_type::iterator, bool> result = counts[lhs].insert(std::make_pair(rule_type::create(rule_type(lhs, rhs)), count));
    if (! result.second)
      result.first->second = cicada::semiring::log(cicada::semiring::traits<weight_type>::exp(result.first->second) + count);
  }
  
  // perform maximum-likelihood estimation...
  // Do we smooth here...?
  
  coarse.clear();
  count_set_type::const_iterator citer_end = counts.end();
  for (count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
    // simple max-like estimation w/o smoothing...
    weight_type sum;
    grammar_type::const_iterator riter_end = citer->second.end();
    for (grammar_type::const_iterator riter = citer->second.begin(); riter != riter_end; ++ riter)
      sum += cicada::semiring::traits<weight_type>::exp(riter->second);
    
    for (grammar_type::const_iterator riter = citer->second.begin(); riter != riter_end; ++ riter)
      coarse[riter->first] = cicada::semiring::log(cicada::semiring::traits<weight_type>::exp(riter->second) / sum);
  }
}

void grammar_counts(const grammar_type& grammar, const lexicon_type& lexicon, expected_counts_type& counts)
{
  count_set_type indexed;
  grammar_type::const_iterator giter_end = grammar.end();
  for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
    indexed[giter->first->lhs].insert(*giter);

  counts.clear();
  counts[goal] = cicada::semiring::traits<weight_type>::one();
  
  expected_counts_type counts_next;
  
  for (int iter = 0; iter < max_iteration; ++ iter) {
    if (debug)
      std::cerr << "iteration: " << (iter + 1) << std::endl;

    counts_next.clear();
    
    expected_counts_type::const_iterator citer_end = counts.end();
    for (expected_counts_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
      count_set_type::const_iterator iiter = indexed.find(citer->first);
      if (iiter == indexed.end())
	continue;
      
      const grammar_type& rules = iiter->second;
      
      grammar_type::const_iterator riter_end = rules.end();
      for (grammar_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	const weight_type weight = cicada::semiring::traits<weight_type>::exp(riter->second) * citer->second;
	
	symbol_set_type::const_iterator siter_end = riter->first->rhs.end();
	for (symbol_set_type::const_iterator siter = riter->first->rhs.begin(); siter != siter_end; ++ siter)
	  if (siter->is_non_terminal())
	    counts_next[*siter] += weight;
      }
    }
    
    // final insertion...!
    for (expected_counts_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
      counts_next.insert(*citer);
    
    counts.swap(counts_next);
    counts_next.clear();
  }
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

void write_grammar(const path_type& prefix, const int order, const grammar_type& grammar)
{  
  typedef std::vector<const grammar_type::value_type*, std::allocator<const grammar_type::value_type*> > sorted_type;
  
  if (grammar.empty()) return;
  
  bool has_suffix_gz  = false;
  bool has_suffix_bz2 = false;
  
  path_type path = prefix;

  if (prefix.extension() == ".gz") {
    path = prefix.parent_path() / prefix.stem();
    has_suffix_gz = true;
  } else if (prefix.extension() == ".bz2") {
    path = prefix.parent_path() / prefix.stem();
    has_suffix_bz2 = true;
  }
  
  path = path.string() + "." + utils::lexical_cast<std::string>(order);
  if (has_suffix_gz)
    path = path.string() + ".gz";
  else if (has_suffix_bz2)
    path = path.string() + ".bz2";

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
    ("output",  po::value<path_type>(&output_file),      "output file (will be augmented by order)")
    
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

