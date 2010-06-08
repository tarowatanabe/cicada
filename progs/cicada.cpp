
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"
#include "cicada/hypergraph.hpp"

#include "cicada/kbest.hpp"

#include "cicada/model.hpp"
#include "cicada/grammar.hpp"

#include "cicada/apply.hpp"
#include "cicada/compose.hpp"
#include "cicada/inside_outside.hpp"
#include "cicada/intersect.hpp"
#include "cicada/binarize.hpp"
#include "cicada/permute.hpp"
#include "cicada/sort.hpp"
#include "cicada/prune.hpp"
#include "cicada/parameter.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/semiring.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>

typedef boost::filesystem::path path_type;

typedef std::vector<std::string, std::allocator<std::string> > grammar_file_set_type;
typedef std::vector<std::string, std::allocator<std::string> > feature_parameter_set_type;

typedef cicada::Symbol          symbol_type;
typedef cicada::Vocab           vocab_type;
typedef cicada::Sentence        sentence_type;
typedef cicada::Lattice         lattice_type;
typedef cicada::Rule            rule_type;
typedef cicada::HyperGraph      hypergraph_type;
typedef cicada::Grammar         grammar_type;
typedef cicada::Model           model_type;
typedef cicada::FeatureFunction feature_function_type;
typedef cicada::Parameter       parameter_type;

typedef cicada::WeightVector<double> weight_set_type;

typedef std::string operation_type;
typedef std::vector<operation_type, std::allocator<operation_type> > operation_set_type;

typedef boost::shared_ptr<hypergraph_type> hypergraph_ptr_type;
typedef std::pair<int, hypergraph_ptr_type> id_hypergraph_type;
typedef std::deque<id_hypergraph_type, std::allocator<id_hypergraph_type> > stack_type;

struct kbest_function
{
  typedef rule_type::feature_set_type feature_set_type;

  typedef cicada::semiring::Logprob<double> value_type;

  kbest_function(const weight_set_type& __weights)
    : weights(__weights) {}

  const weight_set_type& weights;
  
  template <typename Edge>
  value_type operator()(const Edge& edge) const
  {
    return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
  }

  value_type operator()(const feature_set_type& features) const
  {
    return cicada::semiring::traits<value_type>::log(features.dot(weights));
  }

};

struct kbest_function_one
{
  typedef rule_type::feature_set_type feature_set_type;

  typedef cicada::semiring::Logprob<double> value_type;
  
  kbest_function_one(const weight_set_type& __weights) {}

  template <typename Edge>
  value_type operator()(const Edge& edge) const
  {
    return cicada::semiring::traits<value_type>::log(edge.features.dot());
  }

  value_type operator()(const feature_set_type& features) const
  {
    return cicada::semiring::traits<value_type>::log(features.dot());
  }
};


struct kbest_traversal
{
  typedef rule_type::feature_set_type feature_set_type;
  
  typedef boost::tuple<sentence_type, feature_set_type> value_type;
  
  template <typename Edge, typename Iterator>
  void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
  {
    // extract target-yield, features

    boost::get<0>(yield).clear();
    boost::get<1>(yield) = edge.features;
    
    rule_type::symbol_set_type::const_iterator titer_end = edge.rule->target.end();
    for (rule_type::symbol_set_type::const_iterator titer = edge.rule->target.begin(); titer != titer_end; ++ titer)
      if (titer->is_non_terminal()) {
	const int pos = titer->non_terminal_index() - 1;
	boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*(*(first + pos))).begin(), boost::get<0>(*(*(first + pos))).end());
      } else if (*titer != vocab_type::EPSILON)
	boost::get<0>(yield).push_back(*titer);
    
    // collect features...
    for (/**/; first != last; ++ first)
      boost::get<1>(yield) += boost::get<1>(*(*first));
  }
};

struct kbest_filter
{
  kbest_filter(const hypergraph_type& graph) {}
  
  template <typename Node, typename Yield>
  bool operator()(const Node& node, const Yield& yield) const
  {
    return false;
  }
};

struct kbest_filter_unique
{
 #ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<sentence_type, boost::hash<sentence_type>, std::equal_to<sentence_type>, std::allocator<sentence_type> > unique_type;
#else
  typedef sgi::hash_set<sentence_type, boost::hash<sentence_type>, std::equal_to<sentence_type>, std::allocator<sentence_type> > unique_type;
#endif
  typedef std::vector<unique_type, std::allocator<unique_type> > unique_set_type;
 

  kbest_filter_unique(const hypergraph_type& graph) : uniques(graph.nodes.size()) {}
  
  template <typename Node, typename Yield>
  bool operator()(const Node& node, const Yield& yield) const
  {
    unique_set_type& sents = const_cast<unique_set_type&>(uniques);
    unique_type::iterator iter = sents[node.id].find(boost::get<0>(yield));
    if (iter == sents[node.id].end()) {
      sents[node.id].insert(boost::get<0>(yield));
      return false;
    } else
      return true;
  }

  unique_set_type uniques;
};

template <typename Traversal, typename Function, typename Filter>
void kbest_derivations(std::ostream& os,
		       const int id,
		       const hypergraph_type& graph,
		       const int kbest_size,
		       const Traversal& traversal, 
		       const Function& function,
		       const Filter& filter)
{
  cicada::KBest<Traversal, Function, Filter> derivations(graph, kbest_size, traversal, function, filter);
  
  typename Traversal::value_type derivation;
  
  for (int k = 0; k < kbest_size; ++ k) {
    if (! derivations(k, derivation))
      break;

    if (id >= 0)
      os << id << " ||| " << boost::get<0>(derivation) << " |||";
    else
      os << boost::get<0>(derivation) << " |||";
    
    rule_type::feature_set_type::const_iterator fiter_end = boost::get<1>(derivation).end();
    for (rule_type::feature_set_type::const_iterator fiter = boost::get<1>(derivation).begin(); fiter != fiter_end; ++ fiter)
      os << ' ' << fiter->first << '=' << fiter->second;
    os << " ||| ";
    os << function(boost::get<1>(derivation));
    os << '\n';
  }
}

inline
bool true_false(const std::string& token)
{
  if (strcasecmp(token.c_str(), "true") == 0)
    return true;
  if (strcasecmp(token.c_str(), "yes") == 0)
    return true;
  if (atoi(token.c_str()) > 0)
    return true;
  return false;
}

operation_set_type operations;
int debug = 0;

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    stack_type stack;
    
    for (operation_set_type::const_iterator iter = operations.begin(); iter != operations.end(); ++ iter) {
      parameter_type param(*iter);
      
      if (param.name() == "input") {
	bool input_id_mode = (param.find("id") != param.end() ? true_false(param.find("id")->second) : true);
	const std::string path = (param.find("file") != param.end() ? param.find("file")->second : std::string("-"));

	utils::compress_istream is(path, 1024 * 1024);
	
	if (input_id_mode) {
	  int id;
	  std::string sep;
	  hypergraph_ptr_type graph(new hypergraph_type());
	  
	  is >> id >> sep >> *graph;
	  
	  if (! is)
	    throw std::runtime_error("no graph");
	  if (sep != "|||")
	    throw std::runtime_error("hypergraph format error");
	  
	  stack.push_back(std::make_pair(id, graph));
	} else {
	  hypergraph_ptr_type graph(new hypergraph_type());
	  is >> *graph;
	  if (! is)
	    throw std::runtime_error("no graph");
	  
	  stack.push_back(std::make_pair(-1, graph));
	}

      } else if (param.name() == "output") {
	const std::string path = (param.find("file") != param.end() ? param.find("file")->second : std::string("-"));
	
	if (stack.empty())
	  throw std::runtime_error("no graph in the stack!");

	const id_hypergraph_type id_graph = stack.back();
	
	utils::compress_ostream os(path, 1024 * 1024);

	if (id_graph.first >= 0)
	  os << id_graph.first << " ||| " << *id_graph.second << std::endl;
	else
	  os << *id_graph.second << std::endl;
	
      } else if (param.name() == "kbest") {
	const int kbest_size = (param.find("size") != param.end() ? boost::lexical_cast<int>(param.find("size")->second) : 1);
	const std::string path = (param.find("file") != param.end() ? param.find("file")->second : std::string("-"));
	const std::string path_weights = (param.find("weights") != param.end() ? param.find("weights")->second : std::string());
	const bool unique = (param.find("unique") != param.end() ? true_false(param.find("unique")->second) : false);
	
	if (stack.empty())
	  throw std::runtime_error("no graph in the stack!");
	
	const id_hypergraph_type id_graph = stack.back();
	
	weight_set_type weights;

	utils::compress_ostream os(path, 1024 * 1024);

	if (path_weights.empty()) {

	  if (unique)
	    kbest_derivations(os, id_graph.first, *id_graph.second, kbest_size, kbest_traversal(), kbest_function_one(weights), kbest_filter_unique(*id_graph.second));
	  else
	    kbest_derivations(os, id_graph.first, *id_graph.second, kbest_size, kbest_traversal(), kbest_function_one(weights), kbest_filter(*id_graph.second));

	} else {
	  utils::compress_istream is(path_weights);
	  is >> weights;
	  
	  if (unique)
	    kbest_derivations(os, id_graph.first, *id_graph.second, kbest_size, kbest_traversal(), kbest_function(weights), kbest_filter_unique(*id_graph.second));
	  else
	    kbest_derivations(os, id_graph.first, *id_graph.second, kbest_size, kbest_traversal(), kbest_function(weights), kbest_filter(*id_graph.second));
	}
	
	os << std::flush;
      } else if (param.name() == "beam-prune") {
	const double threshold = (param.find("threshold") != param.end() ? boost::lexical_cast<double>(param.find("threshold")->second) : 0.0);
	const std::string path_weights = (param.find("weights") != param.end() ? param.find("weights")->second : std::string());
	
	if (threshold <= 0.0 || 1.0 <= threshold)
	  throw std::runtime_error("pruning is not in a valid range of 0.0 < threshold < 1.0");
	
	if (stack.empty())
	  throw std::runtime_error("no hypergraph in our stack");
	
	weight_set_type weights;
	
	if (path_weights.empty()) {
	  weights.allocate();
	  std::fill(weights.begin(), weights.end(), 1.0);
	} else {
	  utils::compress_istream is(path_weights);
	  is >> weights;
	}
	
	beam_prune(*stack.back().second, weights, 1.0, threshold);
      } else if (param.name() == "unite") {
      
	if (stack.size() < 2)
	  throw std::runtime_error("we need at least two graphs in our stack");
	
	id_hypergraph_type second = stack.back();
	stack.pop_back();
	
	id_hypergraph_type first = stack.back();
	stack.pop_back();
	
	first.second->unite(*second.second);
	
	stack.push_back(first);
      } else if (param.name() == "pop") {
	if (stack.empty())
	  throw std::runtime_error("no graph?");
	
	stack.pop_back();
      }
    } 
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description opts_hidden;
  opts_hidden.add_options()
    ("operation", po::value<operation_set_type>(&operations), "operations");
  
  po::positional_options_description opts_pos;
  opts_pos.add("operation", -1);
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_command;
  
  desc_command.add(opts_command).add(opts_hidden);
  
  po::variables_map variables;
  
  po::store(po::command_line_parser(argc, argv).options(desc_command).positional(opts_pos).run(), variables);
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options] [operations]\n"
	      << opts_command << std::endl;
    exit(0);
  }
}
