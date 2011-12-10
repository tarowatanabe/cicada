//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/generate.hpp>

#include <cicada/operation/normalize.hpp>

#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    
    Normalize::Normalize(const std::string& parameter,
			 const int __debug)
      : base_type("normalize"),
	debug(__debug)
    {
      typedef cicada::Parameter param_type;

      param_type param(parameter);
      if (utils::ipiece(param.name()) != "normalize")
	throw std::runtime_error("this is not a feature normalizer");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "prefix")
	  feature_prefix = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for feature normalizer: " << piter->first << "=" << piter->second << std::endl;
      }
    }
    
    void Normalize::operator()(data_type& data) const
    {
      typedef hypergraph_type::feature_set_type feature_set_type;

      if (! data.hypergraph.is_valid()) return;

      hypergraph_type& hypergraph = data.hypergraph;

      if (debug)
	std::cerr << name << ": " << data.id << std::endl;

      utils::resource start;
	
      hypergraph_type::node_set_type::const_iterator niter_end = hypergraph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = hypergraph.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;

	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  hypergraph_type::edge_type& edge = hypergraph.edges[*eiter];
	    
	  if (! feature_prefix.empty()) {
	    feature_set_type::data_type sum(0.0);
	      
	    feature_set_type::iterator fiter_end = edge.features.end();
	    for (feature_set_type::iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	      if (equal_prefix(feature_prefix, fiter->first))
		sum += fiter->second;
	      
	    if (sum != 0.0)
	      for (feature_set_type::iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
		if (equal_prefix(feature_prefix, fiter->first))
		  fiter->second /= sum;
	  } else {
	    feature_set_type::data_type sum(0.0);
	      
	    feature_set_type::iterator fiter_end = edge.features.end();
	    for (feature_set_type::iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	      sum += fiter->second;
	      
	    if (sum != 0.0)
	      for (feature_set_type::iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
		fiter->second /= sum;
	  }
	}
      }
	
      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << " thread time: " << (end.thread_time() - start.thread_time())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += hypergraph.nodes.size();
      stat.edge += hypergraph.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
    }
    
  };
};
