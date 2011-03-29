//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// learn coarse grammar
// we assume that non-terminals are annotated by @id, and
// simply erasing bits will uncover reduced grammar
//

#include <stdexcept>
#include <vector>
#include <deque>

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

typedef cicada::Rule   rule_type;
typedef cicada::Symbol symbol_type;

path_type input_grammar_file = "-";
path_type input_lexicon_file = "-";
path_type output_prefix;

int max_order = 6;
int max_iteration = 25;

// naive variational bayes for smoothing... otherwise, dirichlet prior
bool variational_bayes_mode = false;
bool maximum_mode = false;

double prior = 0.01;

int threads = 1;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  
  desc.add_options()
    ("input-grammar",   po::value<path_type>(&input_grammar_file),   "input grammar")
    ("input-lexicon",   po::value<path_type>(&input_lexicon_file),   "input lexical rules")
    ("output-prefix",   po::value<path_type>(&output_prefix),        "output prefix")
    
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

