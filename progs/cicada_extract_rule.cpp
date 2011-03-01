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
#include "utils/mathop.hpp"

typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Rule       rule_type;

typedef rule_type::symbol_type     symbol_type;
typedef rule_type::symbol_set_type symbol_set_type;

typedef std::vector<symbol_type, std::allocator<symbol_type> > context_type;

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

double prior = 0.01;
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
      
      symbol_count_type     lhs_counts;
      symbol_set_count_type rhs_counts;
      lhs_rhs_count_type    lhs_rhs_counts;
      rhs_lhs_count_type    rhs_lhs_counts;

      rhs_lhs_count_type left_right_counts;
      rhs_lhs_count_type right_left_counts;
      symbol_set_count_type left_counts;
      symbol_set_count_type right_counts;
      
      context_type context_left;
      context_type context_right;
      
      count_set_type::const_iterator iter_end = counts.end();
      for (count_set_type::const_iterator iter = counts.begin(); iter != iter_end; ++ iter) {
	const rule_type& rule = iter->first;
	const double& count = iter->second;
	
	lhs_rhs_counts[rule.lhs][rule.rhs] += count;
	rhs_lhs_counts[rule.rhs][rule.lhs] += count;
	lhs_counts[rule.lhs] += count;
	rhs_counts[rule.rhs] += count;

	const symbol_type& head_left = rule.rhs.front();
	const symbol_type& head_right = rule.rhs.back();
	
	context_left.clear();
	context_left.push_back(rule.lhs);
	context_left.insert(context_left.end(), rule.rhs.begin(), rule.rhs.end() - 1);
	
	context_right.clear();
	context_right.push_back(rule.lhs);
	context_right.insert(context_right.end(), rule.rhs.begin() + 1, rule.rhs.end());
	
	const symbol_set_type left(context_left.begin(), context_left.end());
	const symbol_set_type right(context_right.begin(), context_right.end());
	
	left_right_counts[left][head_right] += count;
	right_left_counts[right][head_left] += count;
	
	left_counts[left]   += count;
	right_counts[right] += count;
      }
      
      utils::compress_ostream os(output_file);
      os.precision(20);

      if (variational_bayes_mode) {
	lhs_rhs_count_type::const_iterator liter_end = lhs_rhs_counts.end();
	for (lhs_rhs_count_type::const_iterator liter = lhs_rhs_counts.begin(); liter != liter_end; ++ liter) {
	  const symbol_type& lhs = liter->first;
	
	  const double count_lhs = lhs_counts[lhs];
	  const double observed_lhs = liter->second.size();
	  const double factor_lhs = utils::mathop::digamma(count_lhs + observed_lhs * prior);
	
	  symbol_set_count_type::const_iterator riter_end = liter->second.end();
	  for (symbol_set_count_type::const_iterator riter = liter->second.begin(); riter != riter_end; ++ riter) {
	    const symbol_set_type& rhs = riter->first;
	    const double& count = riter->second;
	    
	    const double count_rhs = rhs_counts[rhs];
	    const double observed_rhs = rhs_lhs_counts[rhs].size();
	    const double factor_rhs = utils::mathop::digamma(count_rhs + observed_rhs * prior);
	    
	    const symbol_type& head_left  = rhs.front();
	    const symbol_type& head_right = rhs.back();
	    
	    context_left.clear();
	    context_left.push_back(lhs);
	    context_left.insert(context_left.end(), rhs.begin(), rhs.end() - 1);
	    
	    context_right.clear();
	    context_right.push_back(lhs);
	    context_right.insert(context_right.end(), rhs.begin() + 1, rhs.end());
	    
	    const symbol_set_type left(context_left.begin(), context_left.end());
	    const symbol_set_type right(context_right.begin(), context_right.end());
	    
	    const double count_left  = left_counts[left];
	    const double count_right = right_counts[right];
	    
	    const double observed_left  = left_right_counts[left].size();
	    const double observed_right = right_left_counts[right].size();
	    
	    const double factor_left  = utils::mathop::digamma(count_left + observed_left * prior);
	    const double factor_right = utils::mathop::digamma(count_right + observed_right * prior);
	  
	    os << lhs << " ||| " << rhs << " ||| ||| "
	       << utils::mathop::digamma(count + prior) - factor_lhs << ' '
	       << utils::mathop::digamma(count + prior) - factor_rhs << ' '
	       << utils::mathop::digamma(count + prior) - factor_left << ' '
	       << utils::mathop::digamma(count + prior) - factor_right << '\n';
	  }
	}
      } else {
	lhs_rhs_count_type::const_iterator liter_end = lhs_rhs_counts.end();
	for (lhs_rhs_count_type::const_iterator liter = lhs_rhs_counts.begin(); liter != liter_end; ++ liter) {
	  const symbol_type& lhs = liter->first;
	
	  const double count_lhs = lhs_counts[lhs];
	  const double observed_lhs = liter->second.size();
	  const double factor_lhs = 1.0 / (count_lhs + observed_lhs * prior);
	
	  symbol_set_count_type::const_iterator riter_end = liter->second.end();
	  for (symbol_set_count_type::const_iterator riter = liter->second.begin(); riter != riter_end; ++ riter) {
	    const symbol_set_type& rhs = riter->first;
	    const double& count = riter->second;
	  
	    const double count_rhs = rhs_counts[rhs];
	    const double observed_rhs = rhs_lhs_counts[rhs].size();
	    const double factor_rhs = 1.0 / (count_rhs + observed_rhs * prior);

	    const symbol_type& head_left  = rhs.front();
	    const symbol_type& head_right = rhs.back();
	    
	    context_left.clear();
	    context_left.push_back(lhs);
	    context_left.insert(context_left.end(), rhs.begin(), rhs.end() - 1);
	    
	    context_right.clear();
	    context_right.push_back(lhs);
	    context_right.insert(context_right.end(), rhs.begin() + 1, rhs.end());
	    
	    const symbol_set_type left(context_left.begin(), context_left.end());
	    const symbol_set_type right(context_right.begin(), context_right.end());
	    
	    const double count_left  = left_counts[left];
	    const double count_right = right_counts[right];
	    
	    const double observed_left  = left_right_counts[left].size();
	    const double observed_right = right_left_counts[right].size();

	    const double factor_left  = 1.0 / (count_left + observed_left * prior);
	    const double factor_right = 1.0 / (count_right + observed_right * prior);
	    
	    os << lhs << " ||| " << rhs << " ||| ||| "
	       << utils::mathop::log((count + prior) * factor_lhs) << ' '
	       << utils::mathop::log((count + prior) * factor_rhs) << ' '
	       << utils::mathop::log((count + prior) * factor_left) << ' '
	       << utils::mathop::log((count + prior) * factor_right) << '\n';
	  }
	}
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
    
    ("score",       po::bool_switch(&score_mode),                    "merge and estimate score(s)")
    ("prior",       po::value<double>(&prior)->default_value(prior), "prior")
    ("variational", po::bool_switch(&variational_bayes_mode),        "variational Bayes estimates")
    
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
