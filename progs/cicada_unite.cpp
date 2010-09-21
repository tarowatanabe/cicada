
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

bool individual = false;

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024);

    hypergraph_type merged;
    hypergraph_type hypergraph;

    cicada::Feature feature_confidence(confidence);
    cicada::Feature feature_count(count);
    
    int rank = 1;
    std::string line;
    while (std::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      if (! hypergraph.assign(iter, end))
	throw std::runtime_error("invalid hypergraph format");

      if (! hypergraph.is_valid()) continue;
      
      if (! feature_confidence.empty()) {
	const double conf = 1.0 / (1.0 + rank);
	
	hypergraph_type::edge_set_type::iterator eiter_end = hypergraph.edges.end();
	for (hypergraph_type::edge_set_type::iterator eiter = hypergraph.edges.begin(); eiter != eiter_end; ++ eiter)
	  eiter->features[feature_confidence] = conf;
      } 
      
      if (! feature_count.empty()) {
	hypergraph_type::edge_set_type::iterator eiter_end = hypergraph.edges.end();
	for (hypergraph_type::edge_set_type::iterator eiter = hypergraph.edges.begin(); eiter != eiter_end; ++ eiter)
	  eiter->features[feature_count] = 1;
      }
      
      if (individual)
	os << hypergraph << '\n';
      else
	merged.unite(hypergraph);
      
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
    ("input",  po::value<path_type>(&input_file)->default_value("-"),   "input hypergraphs")
    ("output", po::value<path_type>(&output_file)->default_value("-"),  "output merged hypergraph")
    
    ("confidence", po::value<std::string>(&confidence),    "add confidence weight")
    ("count",      po::value<std::string>(&count),         "add count weight")
    
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
