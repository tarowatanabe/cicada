
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include "cicada_tool_impl.hpp"

path_type operation_file;
operation_set_type operations;
int debug = 0;

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (operations.empty()) {
      if (operation_file.empty())
	operation_file = "-";

      std::string line;
      utils::compress_istream is(operation_file);
      
      while (std::getline(is, line)) {
	typedef boost::tokenizer<utils::space_separator> tokenizer_type;
	
	tokenizer_type tokenizer(line);
	
	operations.clear();
	operations.insert(operations.end(), tokenizer.begin(), tokenizer.end());
	
	process(operations);
      }
    } else
      process(operations);
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
    ("operation-set", po::value<operation_set_type>(&operations), "operations");
  
  po::positional_options_description opts_pos;
  opts_pos.add("operation-set", -1);
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("operation", po::value<path_type>(&operation_file), "operation file (one-line indicate a seriese of operations)")
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
