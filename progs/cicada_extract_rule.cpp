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
#include "utils/sgi_hash_map.hpp"

typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Rule       rule_type;

typedef rule_type::symbol_type     symbol_type;
typedef rule_type::symbol_set_type symbol_set_type;

typedef cicada::semiring::Logprob<double> weight_type;
typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;

#ifdef HAVE_TR1_UNORDERED_MAP
typedef std::tr1::unordered_map<rule_type, double, boost::hash<rule_type>, std::equal_to<rule_type>,
				std::allocator<std::pair<const rule_type, double> > > count_set_type;
#else
typedef sgi::hash_map<rule_type, double, boost::hash<rule_type>, std::equal_to<rule_type>,
		      std::allocator<std::pair<const rule_type, double> > > count_set_type;
#endif

#ifdef HAVE_TR1_UNORDERED_MAP
typedef std::tr1::unordered_map<symbol_type, double, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				std::allocator<std::pair<const symbol_type, double> > > symbol_count_type;
typedef std::tr1::unordered_map<symbol_set_type, double, boost::hash<symbol_set_type>, std::equal_to<symbol_set_type>,
				std::allocator<std::pair<const symbol_set_type, double> > > symbol_set_count_type;

typedef std::tr1::unordered_map<symbol_type, symbol_set_count_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				std::allocator<std::pair<const symbol_type, symbol_set_count_type> > > lhs_rhs_count_type;
typedef std::tr1::unordered_map<symbol_set_type, symbol_count_type, boost::hash<symbol_set_type>, std::equal_to<symbol_set_type>,
				std::allocator<std::pair<const symbol_set_type, symbol_count_type> > > rhs_lhs_count_type;
#else
typedef sgi::hash_map<symbol_type, double, boost::hash<symbol_type>, std::equal_to<symbol_type>,
		      std::allocator<std::pair<const symbol_type, double> > > symbol_count_type;
typedef sgi::hash_map<symbol_set_type, double, boost::hash<symbol_set_type>, std::equal_to<symbol_set_type>,
		      std::allocator<std::pair<const symbol_set_type, double> > > symbol_set_count_type;

typedef sgi::hash_map<symbol_type, symbol_set_count_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
		      std::allocator<std::pair<const symbol_type, symbol_set_count_type> > > lhs_rhs_count_type;
typedef sgi::hash_map<symbol_set_type, symbol_count_type, boost::hash<symbol_set_type>, std::equal_to<symbol_set_type>,
		      std::allocator<std::pair<const symbol_set_type, symbol_count_type> > > rhs_lhs_count_type;
#endif

typedef boost::filesystem::path path_type;

path_type input_file = "-";
path_type output_file = "-";

double prior_terminal = 0.1;
bool score_mode = false;
bool variational_bayes_mode = false;

int debug = 0;

void options(int argc, char** argv);

template <typename Operation>
void process(std::istream& is, Operation op);

struct Output
{
  Output(std::ostream& __os) : os(__os) {}

  void operator()(const rule_type& rule, const double weight)
  {
    os << rule << " ||| ||| " << weight << '\n';
  }
  
  std::ostream& os;
};

struct Counts
{
  Counts(count_set_type& __counts) : counts(__counts) {}

  void operator()(const rule_type& rule, const double weight)
  {
    counts[rule] += weight;
  }
  
  count_set_type& counts;
};

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    utils::compress_istream is(input_file, 1024 * 1024);
    
    if (score_mode) {
      count_set_type counts;
      
      process(is, Counts(counts));
      
      symbol_count_type     lhs;
      symbol_set_count_type rhs;
      lhs_rhs_count_type lhs_rhs;
      rhs_lhs_count_type rhs_lhs;
      
      count_set_type::const_iterator iter_end = counts.end();
      for (count_set_type::const_iterator iter = counts.begin(); iter != iter_end; ++ iter) {
	const rule_type& rule = iter->first;
	const double& count = iter->second;
	
	lhs_rhs[rule.lhs][rule.rhs] += count;
	rhs_lhs[rule.rhs][rule.lhs] += count;
	lhs[rule.lhs] += count;
	rhs[rule.rhs] += count;
      }
      
      utils::compress_ostream os(output_file);
      os.precision(20);
      
      for (count_set_type::const_iterator iter = counts.begin(); iter != iter_end; ++ iter) {
	const rule_type& rule = iter->first;
	const double& count = iter->second;
	
	os << rule << " ||| ||| " << count << '\n';
	
      }
      
      
    } else {
      utils::compress_ostream os(output_file);
      os.precision(20);
      
      process(is, Output(os));
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Operation>
void process(std::istream& is, Operation op)
{
  hypergraph_type graph;
  
  weight_set_type weights_inside;
  weight_set_type weights_edge;
  
  while (is >> graph) {
    if (! graph.is_valid()) continue;
      
    weights_inside.clear();
    weights_edge.clear();
	
    weights_inside.reserve(graph.nodes.size());
    weights_edge.reserve(graph.edges.size());
	
    weights_inside.resize(graph.nodes.size());
    weights_edge.resize(graph.edges.size());
    
    cicada::inside_outside(graph,
			   weights_inside,
			   weights_edge,
			   cicada::operation::weight_function_one<weight_type>(),
			   cicada::operation::weight_function_one<weight_type>());
    
    const weight_type weight_sum = weights_inside.back();
    
    hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
      const hypergraph_type::node_type& node = *niter;
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	
	op(*edge.rule, weights_edge[*eiter] / weight_sum);
      }
    }
  }
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",     po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output",    po::value<path_type>(&output_file)->default_value(output_file), "output")
    
    ("score",          po::bool_switch(&score_mode),                                      "merge and estimate score(s)")
    ("prior-terminal", po::value<double>(&prior_terminal)->default_value(prior_terminal), "prior for terminal")
    ("variational",    po::bool_switch(&variational_bayes_mode),                          "variational Bayes estimates")
    
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
