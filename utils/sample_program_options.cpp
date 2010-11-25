//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//
// -*- encoding: utf-8 -*-

#include <iostream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <utils/program_options.hpp>

std::string op_string;
int op_int;
double op_double;
bool op_bool;


int getoptions(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    if (getoptions(argc, argv) != 0) 
      return 1;

    std::cout << "string: " << op_string << std::endl
	      << "int: " << op_int << std::endl
	      << "double: " << op_double << std::endl
	      << "bool: " << op_bool << std::endl;
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

int getoptions(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("string", po::value<std::string>(&op_string),  "string")
    ("int",    po::value<int>(&op_int),             "int")
    ("double", po::value<double>(&op_double),       "double")
    ("bool",   utils::true_false_switch(&op_bool),  "bool")
    
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    return 1;
  }
  
  return 0;
}
