//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unistd.h>

#include "cicada_impl.hpp"

#include "utils/program_options.hpp"
#include <boost/program_options.hpp>

typedef std::string op_type;
typedef std::vector<op_type, std::allocator<op_type> > op_set_type;

path_type input_file = "-";

bool input_id_mode = false;
bool input_bitext_mode = false;
bool input_lattice_mode = false;
bool input_forest_mode = false;
bool input_span_mode = false;
bool input_directory_mode = false;

std::string symbol_goal         = vocab_type::S;
std::string symbol_non_terminal = vocab_type::X;
path_type   symbol_fallback_file;

grammar_file_set_type grammar_mutable_files;
grammar_file_set_type grammar_static_files;

bool grammar_glue_straight = false;
bool grammar_glue_inverted = false;
bool grammar_insertion = false;
bool grammar_deletion = false;

grammar_file_set_type tree_grammar_mutable_files;
grammar_file_set_type tree_grammar_static_files;

bool tree_grammar_fallback = false;

feature_parameter_set_type feature_parameters;
bool feature_list = false;

op_set_type ops;
bool op_list = false;


int debug = 0;


// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    

    if (feature_list) {
      std::cout << cicada::FeatureFunction::lists();
      return 0;
    }

    if (op_list) {
      std::cout << operation_set_type::lists();
      return 0;
    }

    // read grammars...
    grammar_type grammar;
    const size_t grammar_static_size  = load_grammar<cicada::GrammarStatic>(grammar, grammar_static_files);
    const size_t grammar_mutable_size = load_grammar<cicada::GrammarMutable>(grammar, grammar_mutable_files);
    
    if (debug)
      std::cerr << "loaded static grammar: " << grammar_static_size << std::endl
		<< "loaded mutable grammar: " << grammar_mutable_size << std::endl;
    
    if (grammar_glue_straight || grammar_glue_inverted) {
      if (! symbol_fallback_file.empty()) {
	if (symbol_fallback_file != "-" && ! boost::filesystem::exists(symbol_fallback_file))
	  throw std::runtime_error("invalid fallback non-terminal file: " + symbol_fallback_file.file_string());
	
	utils::compress_istream is(symbol_fallback_file, 1024 * 1024);
	grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarGlue(symbol_goal,
										    symbol_non_terminal,
										    std::istream_iterator<std::string>(is),
										    std::istream_iterator<std::string>(),
										    grammar_glue_straight,
										    grammar_glue_inverted)));
	
      } else
	grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarGlue(symbol_goal,
										    symbol_non_terminal,
										    grammar_glue_straight,
										    grammar_glue_inverted)));
    }

    if (debug)
      std::cerr << "grammar: " << grammar.size() << std::endl;

    tree_grammar_type tree_grammar;
    const size_t tree_grammar_static_size  = load_grammar<cicada::TreeGrammarStatic>(tree_grammar, tree_grammar_static_files);
    const size_t tree_grammar_mutable_size = load_grammar<cicada::TreeGrammarMutable>(tree_grammar, tree_grammar_mutable_files);
    
    if (debug)
      std::cerr << "loaded static tree grammar: " << tree_grammar_static_size << std::endl
		<< "loaded mutable tree grammar: " << tree_grammar_mutable_size << std::endl;
    
    if (debug)
      std::cerr << "tree grammar: " << tree_grammar.size() << std::endl;
    
    // read features...
    model_type model;
    for (feature_parameter_set_type::const_iterator piter = feature_parameters.begin(); piter != feature_parameters.end(); ++ piter)
      model.push_back(feature_function_type::create(*piter));
    model.initialize();

    operation_set_type operations(ops.begin(), ops.end(),
				  model,
				  grammar,
				  tree_grammar,
				  symbol_goal,
				  symbol_non_terminal,
				  grammar_insertion,
				  grammar_deletion,
				  tree_grammar_fallback,
				  input_id_mode || input_directory_mode,
				  input_lattice_mode,
				  input_forest_mode,
				  input_span_mode,
				  input_bitext_mode,
				  false,
				  debug);

    if (! operations.get_output_data().directory.empty()) {
      const path_type& directory = operations.get_output_data().directory;
      
      if (boost::filesystem::exists(directory) && ! boost::filesystem::is_directory(directory))
	boost::filesystem::remove_all(directory);
      
      boost::filesystem::create_directories(directory);
      
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(directory); iter != iter_end; ++ iter)
	boost::filesystem::remove_all(*iter);

      ::sync();
    }
    
    // we will force non directory-input-mode....
    if (input_directory_mode) {
      std::string line;
      
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(input_file); iter != iter_end; ++ iter) {
	utils::compress_istream is(*iter, 1024 * 1024);
	
	if (std::getline(is, line))
	  operations(line);
      }
      
    } else {
      utils::compress_istream is(input_file, 1024 * 1024);
      
      std::string line;
      while (std::getline(is, line))
	operations(line);
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
    
    // options for input/output format
    ("input-id",         po::bool_switch(&input_id_mode),         "id-prefixed input")
    ("input-bitext",     po::bool_switch(&input_bitext_mode),     "target sentence prefixed input")
    ("input-lattice",    po::bool_switch(&input_lattice_mode),    "lattice input")
    ("input-forest",     po::bool_switch(&input_forest_mode),     "forest input")
    ("input-span",       po::bool_switch(&input_span_mode),       "span input")
    ("input-directory",  po::bool_switch(&input_directory_mode),  "input in directory")
    
    // grammar
    ("goal",           po::value<std::string>(&symbol_goal)->default_value(symbol_goal),                 "goal symbol")
    ("non-terminal",   po::value<std::string>(&symbol_non_terminal)->default_value(symbol_non_terminal), "default non-terminal symbol")
    ("fallback",       po::value<path_type>(&symbol_fallback_file),                                      "fallback non-terminal list")
    ("grammar",        po::value<grammar_file_set_type >(&grammar_mutable_files)->composing(),           "grammar file(s)")
    ("grammar-static", po::value<grammar_file_set_type >(&grammar_static_files)->composing(),            "static binary grammar file(s)")
        
    // special handling
    ("grammar-glue-straight", po::bool_switch(&grammar_glue_straight), "add straight hiero glue rule")
    ("grammar-glue-inverted", po::bool_switch(&grammar_glue_inverted), "add inverted hiero glue rule")
    ("grammar-insertion",     po::bool_switch(&grammar_insertion),     "source-to-target transfer grammar")
    ("grammar-deletion",      po::bool_switch(&grammar_deletion),      "source-to-<epsilon> transfer grammar")
    
    // tree-grammar
    ("tree-grammar",          po::value<grammar_file_set_type >(&tree_grammar_mutable_files)->composing(), "tree grammar file(s)")
    ("tree-grammar-static",   po::value<grammar_file_set_type >(&tree_grammar_static_files)->composing(),  "static binary tree grammar file(s)")
    
    // special handling
    ("tree-grammar-fallback", po::bool_switch(&tree_grammar_fallback),                                     "source-to-target transfer tree grammar")
    
    // models...
    ("feature-function",      po::value<feature_parameter_set_type >(&feature_parameters)->composing(), "feature function(s)")
    ("feature-function-list", po::bool_switch(&feature_list),                                           "list of available feature function(s)")

    //operatins...
    ("operation",      po::value<op_set_type>(&ops)->composing(), "operations")
    ("operation-list", po::bool_switch(&op_list),                 "list of available operation(s)");

  

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
    const path_type path_config = variables["config"].as<path_type>();
    if (! boost::filesystem::exists(path_config))
      throw std::runtime_error("no config file: " + path_config.file_string());
    
    utils::compress_istream is(path_config);
    po::store(po::parse_config_file(is, desc_config), variables);
  }
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc_command << std::endl;
    exit(0);
  }
}
