
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

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/semiring.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"

#include "cicada_impl.hpp"

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

typedef cicada::WeightVector<double> weight_set_type;

struct weight_set_function
{
  typedef cicada::semiring::Logprob<double> value_type;

  weight_set_function(const weight_set_type& __weights)
    : weights(__weights) {}

  const weight_set_type& weights;
  
  template <typename FeatureSet>
  value_type operator()(const FeatureSet& x) const
  {
    return cicada::semiring::traits<value_type>::log(x.dot(weights));
  }
};


struct weight_set_function_one
{
  typedef cicada::semiring::Logprob<double> value_type;

  weight_set_function_one(const weight_set_type& __weights) {}
  
  template <typename FeatureSet>
  value_type operator()(const FeatureSet& x) const
  {
    return cicada::semiring::traits<value_type>::log(x.dot());
  }
};


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
		       const size_t id,
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
    
    os << id << " ||| " << boost::get<0>(derivation) << " |||";
    rule_type::feature_set_type::const_iterator fiter_end = boost::get<1>(derivation).end();
    for (rule_type::feature_set_type::const_iterator fiter = boost::get<1>(derivation).begin(); fiter != fiter_end; ++ fiter)
      os << ' ' << fiter->first << '=' << fiter->second;
    os << " ||| ";
    os << function(boost::get<1>(derivation));
    os << '\n';
  }
}

path_type input_file = "-";
path_type output_file = "-";

bool input_id_mode = false;
bool input_bitext_mode = false;
bool input_lattice_mode = false;
bool input_forest_mode = false;
bool input_directory_mode = false;
bool output_forest_mode = false;
bool output_directory_mode = false;

int kbest_size = 1;
bool kbest_unique = false;

std::string symbol_goal         = vocab_type::S;
std::string symbol_non_terminal = vocab_type::X;

grammar_file_set_type grammar_mutable_files;
grammar_file_set_type grammar_static_files;

bool grammar_glue_straight = false;
bool grammar_glue_inverted = false;
bool grammar_insertion = false;
bool grammar_deletion = false;

feature_parameter_set_type feature_parameters;
path_type                  feature_weights_file;
bool feature_weights_one = false;
bool feature_list = false;

bool binarize_graph = false;
int  binarize_size = 2;

bool permute_graph = false;
int  permute_size = 3;

bool apply_cube = false;
bool apply_full = false;
int  cube_size = 200;

double prune_beam = 0.0;

int debug = 0;


// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (input_lattice_mode && input_forest_mode)
      throw std::runtime_error("input can be sentence, lattice or forest");
    
    if (apply_cube && apply_full)
      throw std::runtime_error("intersection can be either cube or full (default dube)");

    if (feature_list) {
      std::cout << cicada::FeatureFunction::lists();
      return 0;
    }

    // read grammars...
    grammar_type grammar;
    for (grammar_file_set_type::const_iterator fiter = grammar_static_files.begin(); fiter != grammar_static_files.end(); ++ fiter)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarStatic(*fiter)));

    if (debug)
      std::cerr << "loaded mutable grammar: " << grammar.size() << std::endl;
    
    for (grammar_file_set_type::const_iterator fiter = grammar_mutable_files.begin(); fiter != grammar_mutable_files.end(); ++ fiter)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarMutable(*fiter)));

    if (debug)
      std::cerr << "loaded static grammar: " << grammar.size() << std::endl;
    
    if (grammar_glue_straight || grammar_glue_inverted)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarGlue(symbol_goal,
										  symbol_non_terminal,
										  grammar_glue_straight,
										  grammar_glue_inverted)));

    if (debug)
      std::cerr << "grammar: " << grammar.size() << std::endl;
    
    // read features...
    model_type model;
    for (feature_parameter_set_type::const_iterator piter = feature_parameters.begin(); piter != feature_parameters.end(); ++ piter)
      model.push_back(feature_function_type::create(*piter));

    if (debug)
      std::cerr << "feature functions: " << model.size() << std::endl;

    // read parameters...
    weight_set_type weights;
    if (! feature_weights_file.empty()) {
      if (feature_weights_file != "-" && ! boost::filesystem::exists(feature_weights_file))
	throw std::runtime_error("no feture weights?" + feature_weights_file.file_string());
      
      if (feature_weights_one)
	throw std::runtime_error("feature weights file supplied but you have enabled one-initialized weights");
      
      utils::compress_istream is(feature_weights_file);
      is >> weights;
    }
    
    // we will force non directory-input-mode....
    input_directory_mode = false;
    
    boost::shared_ptr<utils::compress_istream> is(! input_directory_mode  ? new utils::compress_istream(input_file) : 0);
    boost::shared_ptr<utils::compress_ostream> os(! output_directory_mode ? new utils::compress_ostream(output_file, 1024 * 1024) : 0);

    if (output_directory_mode) {
      
      if (boost::filesystem::exists(output_file) && ! boost::filesystem::is_directory(output_file))
	boost::filesystem::remove_all(output_file);
      
      boost::filesystem::create_directories(output_file);
      
      // remove all the files..
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(output_file); iter != iter_end; ++ iter)
	boost::filesystem::remove_all(*iter);
    }
    
    std::string     line;
    sentence_type   target_sentence;
    lattice_type    target;
    sentence_type   sentence;
    lattice_type    lattice;
    hypergraph_type hypergraph;
    hypergraph_type hypergraph_composed;
    hypergraph_type hypergraph_applied;
    
    size_t id = 0;
    while (std::getline(*is, line)) {
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      if (input_id_mode)
	if (! parse_id(id, iter, end))
	  throw std::runtime_error("invalid id-prefixed format");
      
      if (input_lattice_mode) {
	if (! lattice.assign(iter, end))
	  throw std::runtime_error("invalid lattive format");
      } else if (input_forest_mode) {
	if (! hypergraph.assign(iter, end))
	  throw std::runtime_error("invalid hypergraph format");
      } else {
	if (! sentence.assign(iter, end))
	  throw std::runtime_error("invalid sentence format");
	
	lattice = lattice_type(sentence);
      }
      
      if (input_bitext_mode) {
	if (! parse_separator(iter, end))
	  throw std::runtime_error("no separator?");
	
	if (! target_sentence.assign(iter, end))
	  throw std::runtime_error("invalid sentence format");
	
	target = lattice_type(target_sentence);
      }
      
      if (iter != end)
	throw std::runtime_error("invalid input format");
      
      if (lattice.empty() && ! hypergraph.is_valid()) {
	++ id;
	continue;
      }

      
      grammar_type grammar_translation(grammar);

      if (input_forest_mode) {
	if (grammar_insertion)
	  grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(hypergraph, symbol_non_terminal)));
	if (grammar_deletion)
	  grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(hypergraph, symbol_non_terminal)));
      } else {
	if (grammar_insertion)
	  grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(lattice, symbol_non_terminal)));
	if (grammar_deletion)
	  grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(lattice, symbol_non_terminal)));
      }

      if (debug)
	if (input_forest_mode)
	  std::cerr << "# of nodes: " << hypergraph.nodes.size()
		    << " # of edges: " << hypergraph.edges.size()
		    << " valid? " << (hypergraph.is_valid() ? "true" : "false")
		    << std::endl;

      if (input_forest_mode && binarize_graph) {
	hypergraph_type hypergraph_binarized;

	if (debug)
	  std::cerr << "binarization" << std::endl;
	
	utils::resource binarize_start;
	
	cicada::binarize(hypergraph, hypergraph_binarized, binarize_size);
	
	utils::resource binarize_end;
	
	if (debug)
	  std::cerr << "binarize cpu time: " << (binarize_end.cpu_time() - binarize_start.cpu_time())
		    << " user time: " << (binarize_end.user_time() - binarize_start.user_time())
		    << std::endl;
	
	if (debug)
	  std::cerr << "# of nodes: " << hypergraph_binarized.nodes.size()
		    << " # of edges: " << hypergraph_binarized.edges.size()
		    << " valid? " << (hypergraph_binarized.is_valid() ? "true" : "false")
		    << std::endl;
	
	hypergraph.swap(hypergraph_binarized);
      }
      
      if (input_forest_mode && permute_graph) {
	
	hypergraph_type hypergraph_permuted;

	if (debug)
	  std::cerr << "permutation" << std::endl;
	
	utils::resource permute_start;
	
	cicada::permute(hypergraph, hypergraph_permuted, permute_size);
	
	utils::resource permute_end;
	
	if (debug)
	  std::cerr << "permute cpu time: " << (permute_end.cpu_time() - permute_start.cpu_time())
		    << " user time: " << (permute_end.user_time() - permute_start.user_time())
		    << std::endl;
	
	if (debug)
	  std::cerr << "# of nodes: " << hypergraph_permuted.nodes.size()
		    << " # of edges: " << hypergraph_permuted.edges.size()
		    << " valid? " << (hypergraph_permuted.is_valid() ? "true" : "false")
		    << std::endl;
	
	hypergraph.swap(hypergraph_permuted);
      }


      if (debug)
	std::cerr << "composition" << std::endl;

      utils::resource compose_start;
      
      if (input_forest_mode) {
	// we assume cfg-fst composition
	
	cicada::compose_earley(grammar_translation, hypergraph, hypergraph_composed);
      } else {
	// we assume synchronous-cfg composition
	
	cicada::compose_cky(symbol_goal, grammar_translation, lattice, hypergraph_composed);
      }

      utils::resource compose_end;
      
      if (debug)
	std::cerr << "compose cpu time: " << (compose_end.cpu_time() - compose_start.cpu_time())
		  << " user time: " << (compose_end.user_time() - compose_start.user_time())
		  << std::endl;

      if (debug)
	std::cerr << "# of nodes: " << hypergraph_composed.nodes.size()
		  << " # of edges: " << hypergraph_composed.edges.size()
		  << " valid? " << (hypergraph_composed.is_valid() ? "true" : "false")
		  << std::endl;
      
      if (! model.empty()) {

	if (debug)
	  std::cerr << "apply features" << std::endl;
	
	utils::resource apply_start;
	
	if (apply_full) {
	  if (feature_weights_one)
	    cicada::apply_exact<weight_set_function_one>(model, hypergraph_composed, hypergraph_applied, weight_set_function_one(weights), cube_size);
	  else
	    cicada::apply_exact<weight_set_function>(model, hypergraph_composed, hypergraph_applied, weight_set_function(weights), cube_size);
	} else {
	  if (feature_weights_one)
	    cicada::apply_cube_prune<weight_set_function_one>(model, hypergraph_composed, hypergraph_applied, weight_set_function_one(weights), cube_size);
	  else
	    cicada::apply_cube_prune<weight_set_function>(model, hypergraph_composed, hypergraph_applied, weight_set_function(weights), cube_size);
	}
	
	utils::resource apply_end;
	
	if (debug)
	  std::cerr << "apply cpu time: " << (apply_end.cpu_time() - apply_start.cpu_time())
		    << " user time: " << (apply_end.user_time() - apply_start.user_time())
		    << std::endl;

	if (debug)
	  std::cerr << "# of nodes: " << hypergraph_applied.nodes.size()
		    << " # of edges: " << hypergraph_applied.edges.size()
		    << " valid? " << (hypergraph_applied.is_valid() ? "true" : "false")
		    << std::endl;
      } else
	hypergraph_applied = hypergraph_composed;

      
      if (0.0 < prune_beam && prune_beam < 1.0) {

	hypergraph_type hypergraph_pruned;
	
	utils::resource prune_start;
	
	beam_prune(hypergraph_applied, hypergraph_pruned, weights, 1.0, prune_beam);
	
	utils::resource prune_end;
	
	if (debug)
	std::cerr << "prune cpu time: " << (prune_end.cpu_time() - prune_start.cpu_time())
		  << " user time: " << (prune_end.user_time() - prune_start.user_time())
		  << std::endl;

	if (debug)
	  std::cerr << "# of nodes: " << hypergraph_pruned.nodes.size()
		    << " # of edges: " << hypergraph_pruned.edges.size()
		    << " valid? " << (hypergraph_pruned.is_valid() ? "true" : "false")
		    << std::endl;

	hypergraph_applied.swap(hypergraph_pruned);
      }
      
      if (input_bitext_mode) {
	
	hypergraph_type hypergraph_intersected;
	
	utils::resource intersect_start;
	
	intersect(hypergraph_applied, target, hypergraph_intersected);
	
	utils::resource intersect_end;
	
	if (debug)
	  std::cerr << "intersect cpu time: " << (intersect_end.cpu_time() - intersect_start.cpu_time())
		    << " user time: " << (intersect_end.user_time() - intersect_start.user_time())
		    << std::endl;
	
	if (debug)
	  std::cerr << "# of nodes: " << hypergraph_intersected.nodes.size()
		    << " # of edges: " << hypergraph_intersected.edges.size()
		    << " valid? " << (hypergraph_intersected.is_valid() ? "true" : "false")
		    << std::endl;
	
	hypergraph_applied.swap(hypergraph_intersected);
      }
      
      if (hypergraph_applied.is_valid()) {
	
	if (output_directory_mode) {
	  const path_type path = path_type(output_file) / (boost::lexical_cast<std::string>(id) + ".gz");
	  
	  os.reset(new utils::compress_ostream(path, 1024 * 1024));
	}
	
	if (output_forest_mode)
	  *os << id << " ||| " << hypergraph_applied << '\n';
	else {
	  // extract k-best ...
	  
	  if (feature_weights_one) {
	    if (kbest_unique)
	      kbest_derivations(*os, id, hypergraph_applied, kbest_size, kbest_traversal(), kbest_function_one(weights), kbest_filter_unique(hypergraph_applied));
	    else
	      kbest_derivations(*os, id, hypergraph_applied, kbest_size, kbest_traversal(), kbest_function_one(weights), kbest_filter(hypergraph_applied));
	  } else {
	    if (kbest_unique)
	      kbest_derivations(*os, id, hypergraph_applied, kbest_size, kbest_traversal(), kbest_function(weights), kbest_filter_unique(hypergraph_applied));
	    else
	      kbest_derivations(*os, id, hypergraph_applied, kbest_size, kbest_traversal(), kbest_function(weights), kbest_filter(hypergraph_applied));
	  }
	}
      }
      
      ++ id;

      *os << std::flush;
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

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    // options for input/output format
    ("input-id",         po::bool_switch(&input_id_mode),         "id-prefixed input")
    ("input-bitext",     po::bool_switch(&input_bitext_mode),     "target sentence prefixed input")
    ("input-lattice",    po::bool_switch(&input_lattice_mode),    "lattice input")
    ("input-forest",     po::bool_switch(&input_forest_mode),     "forest input")
    ("input-directory",  po::bool_switch(&input_directory_mode),  "input in directory")
    ("output-forest",    po::bool_switch(&output_forest_mode),    "forest output")
    ("output-directory", po::bool_switch(&output_directory_mode), "output in directory")
    
    // k-best derivation output
    ("kbest",        po::value<int>(&kbest_size),    "k-best derivation")
    ("kbest-unique", po::bool_switch(&kbest_unique), "unique k-best")
    
    // grammar
    ("goal",           po::value<std::string>(&symbol_goal)->default_value(symbol_goal),                 "goal symbol")
    ("non-terminal",   po::value<std::string>(&symbol_non_terminal)->default_value(symbol_non_terminal), "default non-terminal symbol")
    ("grammar",        po::value<grammar_file_set_type >(&grammar_mutable_files),                        "grammar file(s)")
    ("grammar-static", po::value<grammar_file_set_type >(&grammar_static_files),                         "static binary grammar file(s)")
    
    // special handling
    ("grammar-glue-straight", po::bool_switch(&grammar_glue_straight), "add straight hiero glue rule")
    ("grammar-glue-inverted", po::bool_switch(&grammar_glue_inverted), "add inverted hiero glue rule")
    ("grammar-insertion",     po::bool_switch(&grammar_insertion),     "source-to-target transfer grammar")
    ("grammar-deletion",      po::bool_switch(&grammar_deletion),      "source-to-<epsilon> transfer grammar")
    
    // models...
    ("feature-function",      po::value<feature_parameter_set_type >(&feature_parameters), "feature function(s)")
    ("feature-weights",       po::value<path_type>(&feature_weights_file),                 "feature weights")
    ("feature-weights-one",   po::bool_switch(&feature_weights_one),                       "one initialized feature weights")
    ("feature-function-list", po::bool_switch(&feature_list),                              "list of available feature function(s)")
    
    // binarization
    ("binarize",      po::bool_switch(&binarize_graph),                             "perform hypergraph binarization")
    ("binarize-size", po::value<int>(&binarize_size)->default_value(binarize_size), "binarization size (<=2 for all binzarization)")
    
    // permutation...
    ("permute",      po::bool_switch(&permute_graph),                            "perform hypergraph permutation")
    ("permute-size", po::value<int>(&permute_size)->default_value(permute_size), "permutation size (zero for no-permutation, negative for all permutation)")
    
    // application strategy
    ("apply-cube", po::bool_switch(&apply_cube),                         "feature application by cube-pruning")
    ("apply-full", po::bool_switch(&apply_full),                         "full feature application")
    ("cube-size",  po::value<int>(&cube_size)->default_value(cube_size), "cube-size for cube prunning")
    
    // beam pruning
    ("prune-beam", po::value<double>(&prune_beam),  "beam pruning (0.0 < threshold < 1.0)")
    
    ;

  

  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("config", po::value<path_type>(), "configuration file")
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_config;
  po::options_description desc_command;
  
  desc_config.add(opts_config);
  desc_command.add(opts_config).add(opts_command);
  
  po::variables_map variables;

  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  if (variables.count("config")) {
    utils::compress_istream is(variables["config"].as<path_type>());
    po::store(po::parse_config_file(is, desc_config), variables);
  }
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc_command << std::endl;
    exit(0);
  }
}
