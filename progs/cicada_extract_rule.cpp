//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>

#include <cicada/hypergraph.hpp>
#include <cicada/inside_outside.hpp>
#include <cicada/operation/functional.hpp>
#include <cicada/semiring.hpp>

#include <boost/filesystem.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"

typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Rule       rule_type;

typedef cicada::semiring::Logprob<double> weight_type;
typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;

typedef boost::filesystem::path path_type;

path_type input_file = "-";
path_type output_file = "-";

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    hypergraph_type graph;

    weight_set_type weights_inside;
    weight_set_type weights_edge;
  
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file);
    os.precision(20);
    
    while (is >> graph) {
      if (! graph.is_valid()) continue;
      
      weights_inside.clear();
      weights_edge.clear();
      
      weights_inside.reserve(graph.nodes.size());
      weights_edge.reserve(graph.edges.size());
      
      weights_inside.resize(graph.nodes.size());
      weights_edge.resize(graph.edges.size());

      cicada::inside_outside(graph, weights_inside, weights_edge, cicada::operation::weight_function_one<weight_type>(), cicada::operation::weight_function_one<weight_type>());

      const weight_type weight_sum = weights_inside.back();
      
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  os << *edge.rule << " ||| ||| " << static_cast<double>(weights_edge[*eiter] / weight_sum) << '\n';
	}
      }
      
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
  
  po::options_description desc("options");
  desc.add_options()
    ("input",     po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output",    po::value<path_type>(&output_file)->default_value(output_file), "output")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
        
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}
