
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include "cicada_impl.hpp"

#include "utils/program_options.hpp"

#include <boost/program_options.hpp>

path_type input_file = "-";
path_type output_file = "-";

std::string confidence;
std::string count;
double count_weight = 1.0;

bool individual = false;

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024);
    
    lattice_type merged;
    lattice_type lattice;

    cicada::Feature feature_confidence(confidence);
    cicada::Feature feature_count(count);
    
    int rank = 1;
    std::string line;
    while (std::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      if (! lattice.assign(iter, end))
	throw std::runtime_error("invalid hypergraph format");
      
      if (! feature_confidence.empty()) {
	const double conf = 1.0 / (1.0 + rank);
	
	lattice_type::iterator liter_end = lattice.end();
	for (lattice_type::iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	  lattice_type::arc_set_type& arcs = *liter;
	  
	  for (lattice_type::arc_set_type::iterator aiter = arcs.begin(); aiter != arcs.end(); ++ aiter)
	    aiter->features[feature_confidence] = conf;
	}
      } 
      
      if (! feature_count.empty()) {
	lattice_type::iterator liter_end = lattice.end();
	for (lattice_type::iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	  lattice_type::arc_set_type& arcs = *liter;
	  
	  for (lattice_type::arc_set_type::iterator aiter = arcs.begin(); aiter != arcs.end(); ++ aiter)
	    aiter->features[feature_count] = count_weight;
	}
      }
      
      if (individual)
	os << lattice << '\n';
      else {
	// perform merging...
	
      }
	
      
      ++ rank;
    }
    
    if (! individual)
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
    
    ("confidence",   po::value<std::string>(&confidence),    "add confidence weight feature name")
    ("count",        po::value<std::string>(&count),         "add count weight feature name")
    ("count-weight", po::value<double>(&count_weight),       "count weight")
    
    ("individual", po::bool_switch(&individual), "no merging")
    
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
