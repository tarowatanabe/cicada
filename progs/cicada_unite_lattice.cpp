#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <set>

#include "cicada_impl.hpp"

#include "utils/sgi_hash_map.hpp"
#include "utils/program_options.hpp"

#include <boost/program_options.hpp>

#include <google/dense_hash_set>

path_type input_file = "-";
path_type output_file = "-";

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024);
    
    lattice_type merged;
    lattice_type merged_new;
    lattice_type lattice;
    
    std::string line;
    while (std::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      if (! lattice.assign(iter, end))
	throw std::runtime_error("invalid hypergraph format");
      
      if (lattice.empty()) continue;
      
      if (merged.empty())
	merged.swap(lattice);
      else {
	
	// we simply concatenate together...!
	//
	// new-root -- merged -- new-node ------------- new-goal
	//          ------------          -- lattice -- 
	//
	
	merged_new.clear();
	merged_new.push_back(lattice_type::arc_set_type());
	merged_new.back().push_back(lattice_type::arc_type(vocab_type::EPSILON, feature_set_type(), 1));
	merged_new.back().push_back(lattice_type::arc_type(vocab_type::EPSILON, feature_set_type(), merged.size() + 2));
	
	for (int i = 0; i != merged.size(); ++ i)
	  merged_new.push_back(merged[i]);
	
	merged_new.push_back(lattice_type::arc_set_type());
	merged_new.back().push_back(lattice_type::arc_type(vocab_type::EPSILON, feature_set_type(), lattice.size() + 1));
	
	for (int i = 0; i != lattice.size(); ++ i)
	  merged_new.push_back(lattice[i]);
	
	merged.swap(merged_new);
	merged_new.clear();
      }
    }
    
    os << merged << '\n';
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

  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value("-"),   "input lattices")
    ("output", po::value<path_type>(&output_file)->default_value("-"),  "output merged lattice")
    
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
