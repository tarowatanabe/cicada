
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"

#include "cicada/kbest.hpp"

#include "cicada/model.hpp"
#include "cicada/grammar.hpp"

#include "cicada/apply.hpp"
#include "cicada/compose.hpp"
#include "cicada/inside_outside.hpp"
#include "cicada/intersect.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/semiring.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>

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

struct kbest_function
{
  typedef cicada::semiring::Logprob<double> value_type;

  kbest_function(const weight_set_type& __weights)
    : weights(__weights) {}

  const weight_set_type& weights;
  
  template <typename Edge>
  value_type operator()(const Edge& edge) const
  {
    return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
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
      } else
	boost::get<0>(yield).push_back(*titer);
    
    // collect features...
    for (/**/; first != last; ++ first)
      boost::get<1>(yield) += boost::get<1>(*(*first));
  }
};


path_type input_file = "-";
path_type output_file = "-";

bool input_lattice_mode = false;
bool input_forest_mode = false;
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
bool feature_list = false;

bool intersection_cube = false;
bool intersection_full = false;
int  cube_size = 200;

int debug = 0;


// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (input_lattice_mode && input_forest_mode)
      throw std::runtime_error("input can be sentence, lattice or forest");
    
    if (intersection_cube && intersection_full)
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
      
      utils::compress_istream is(feature_weights_file);
      is >> weights;
    }
    
    
    utils::compress_istream is(input_file);
    utils::compress_ostream os(output_file, 1024 * 1024);

    std::string     line;
    lattice_type    lattice;
    hypergraph_type hypergraph;
    hypergraph_type hypergraph_composed;
    hypergraph_type hypergraph_applied;
    
    size_t id = 0;
    while (1) {
      if (input_lattice_mode)
	is >> lattice;
      else if (input_forest_mode)
	is >> hypergraph;
      else {
	sentence_type sentence;
	is >> sentence;
	if (is)
	  lattice = lattice_type(sentence);
      }
      
      if (! is) break;
      
      grammar_type grammar_translation(grammar);
      if (grammar_insertion)
	grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(lattice, symbol_non_terminal)));
      if (grammar_deletion)
	grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(lattice, symbol_non_terminal)));

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
		  << " valid? " << (hypergraph_composed.goal != hypergraph_type::invalid)
		  << std::endl;
      
      if (debug)
	std::cerr << "apply features" << std::endl;

      utils::resource apply_start;

      cicada::apply_cube_prune<weight_set_function>(model, hypergraph_composed, hypergraph_applied, weight_set_function(weights), cube_size);

      utils::resource apply_end;
      
      if (debug)
	std::cerr << "apply cpu time: " << (apply_end.cpu_time() - apply_start.cpu_time())
		  << " user time: " << (apply_end.user_time() - apply_start.user_time())
		  << std::endl;

      if (debug)
	std::cerr << "# of nodes: " << hypergraph_applied.nodes.size()
		  << " # of edges: " << hypergraph_applied.edges.size()
		  << " valid? " << (hypergraph_applied.goal != hypergraph_type::invalid)
		  << std::endl;
      
      if (output_forest_mode)
	os << hypergraph_applied << '\n';
      else {
	// extract k-best ...
	cicada::KBest<kbest_traversal, kbest_function> kbest_derivations(hypergraph_applied,
									 kbest_size,
									 kbest_traversal(),
									 kbest_function(weights));
	
	kbest_traversal::value_type derivation;
	for (int k = 0; k < kbest_size; ++ k) {
	  if (! kbest_derivations(k, derivation))
	    break;
	  
	  os << id << " ||| " << boost::get<0>(derivation) << " |||";
	  rule_type::feature_set_type::const_iterator fiter_end = boost::get<1>(derivation).end();
	  for (rule_type::feature_set_type::const_iterator fiter = boost::get<1>(derivation).begin(); fiter != fiter_end; ++ fiter)
	    os << ' ' << fiter->first << '=' << fiter->second;
	  os << " ||| ";
	  os << boost::get<1>(derivation).dot(weights);
	  os << '\n';
	}
	++ id;
	
      }
      
      os << std::flush;
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
    ("input-lattice",    po::bool_switch(&input_lattice_mode),    "lattice input")
    ("input-forest",     po::bool_switch(&input_forest_mode),     "forest input")
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
    ("feature-weights",       po::value<path_type>(&feature_weights_file),                  "feature weights")
    ("feature-function-list", po::bool_switch(&feature_list),                              "list of available feature function(s)")
    
    // intersection strategy
    ("intersection-cube", po::bool_switch(&intersection_cube),                  "intersetion by cube-pruning")
    ("intersection-full", po::bool_switch(&intersection_full),                  "full intersection")
    ("cube-size",         po::value<int>(&cube_size)->default_value(cube_size), "cube-size for cube prunning");

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
