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
void read_grammar(const path_type& path, grammar_type& grammar);
void read_lexicon(const path_type& path, lexicon_type& lexicon);

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
    
    for (int order = max_iteration - 1; order >= 0; -- order) {
      
      
    }
    
    // finally reduce to minus-grammar:
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
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
    if (! parse_rule(line, rule, logrpob)) continue;
    
    grammar[rule_type::create(rule)] = logprob;
  }
  
}

void read_lexicon(const path_type& path, lexicon_type& lexicon)
{
  utils::compress_istream is(path, 1024 * 1024);

  utils::compress_istream is(path, 1024 * 1024);
  
  std::string line;
  rule_type rule;
  double logprob;
  
  while (std::getline(is, line)) {
    if (! parse_rule(line, rule, logrpob)) continue;
    
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

