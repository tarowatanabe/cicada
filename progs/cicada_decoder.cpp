
#include <iostream>
#include <vector>
#include <string>

#include "cicada/kbest.hpp"

#include "cicada/model.hpp"
#include "cicada/grammar.hpp"

#include "cicada/apply.hpp"
#include "cicada/compose.hpp"
#include "cicada/inside_outside.hpp"
#include "cicada/intersect.hpp"

#include "cicada/feature_function.hpp"

#include "utils/program_options.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

typedef boost::filesystem::path path_type;

typedef boost::program_options::variables_map variable_set_type;

typedef cicada::Vocab vocab_type;


// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv, variable_set_type& variables);

int main(int argc, char ** argv)
{
  try {
    variable_set_type variables;
    
    options(argc, argv, variables);
    
    if (variables.count("feature-function-list")) {
      std::cout << cicada::FeatureFunction::lists();
      return 0;
    }
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}


void options(int argc, char** argv, variable_set_type& variables)
{
  namespace po = boost::program_options;

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("input",            po::value<path_type>()->default_value("-"), "input file")
    ("output",           po::value<path_type>()->default_value("-"), "output file")
    
    // options for input/output format
    ("input-lattice",    po::bool_switch(), "lattice input")
    ("input-forest",     po::bool_switch(), "forest input")
    ("output-forest",    po::bool_switch(), "forest output")
    ("output-directory", po::bool_switch(), "output in directory")
    
    // k-best derivation output
    ("kbest",        po::value<int>(),  "k-best derivation")
    ("kbest-unique", po::bool_switch(), "unique k-best")
    
    // grammar
    ("goal",                  po::value<std::string>()->default_value(vocab_type::S), "goal symbol")
    ("non-terminal",          po::value<std::string>()->default_value(vocab_type::X), "default non-terminal symbol")
    ("gramamr",               po::value<std::vector<std::string> >(),                 "grammar file(s)")
    ("gramamr-static",        po::value<std::vector<std::string> >(),                 "static binary grammar file(s)")
    
    // special handling
    ("grammar-glue-straight", po::bool_switch(), "add straight hiero glue rule")
    ("grammar-glue-inverted", po::bool_switch(), "add inverted hiero glue rule")
    ("grammar-insertion",     po::bool_switch(), "source-to-target transfer grammar")
    ("grammar-deletion",      po::bool_switch(), "source-to-<epsilon> transfer grammar")
    
    // models...
    ("feature-function",      po::value<std::vector<std::string> >(), "feature function(s)")
    ("feature-function-list", po::bool_switch(),                      "list of available feature function(s)")
    
    // intersection strategy
    ("intersection-cube", po::bool_switch(),                    "intersetion by cube-pruning")
    ("intersection-full", po::bool_switch(),                    "full intersection")
    ("cube-size",         po::value<int>()->default_value(200), "cube-size for cube prunning");

  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("config", po::value<path_type>(), "configuration file")
    ("help", "help message");
  
  po::options_description desc_config;
  po::options_description desc_command;
  
  desc_config.add(opts_config);
  desc_command.add(opts_config).add(opts_command);
  
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
