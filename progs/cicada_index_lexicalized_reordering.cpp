
#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>

#include "cicada/lexicalized_reordering.hpp"

#include "utils/program_options.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>


typedef boost::program_options::variables_map variable_set_type;

void options(int argc, char** argv, variable_set_type& variables);

int main(int argc, char** argv)
{
  try {
    variable_set_type variables;
    options(argc, argv, variables);
    
    if (! variables.count("input"))
      throw std::runtime_error("no input file");
    if (! variables.count("output"))
      throw std::runtime_error("no output file");
    
    const std::string input_path  = variables["input"].as<std::string>();
    const std::string output_path = variables["output"].as<std::string>();
    
    cicada::LexicalizedReordering model(input_path);
    
    if (variables.count("quantize") && variables["quantize"].as<bool>())
      model.quantize();

    model.write(output_path);
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
  
  po::options_description desc("options");
  desc.add_options()
    ("input",  po::value<std::string>()->default_value("-"),   "input in text format")
    ("output", po::value<std::string>(), "output in binary format")

    ("quantize", utils::true_false_switch(), "perform quantization")
    
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
