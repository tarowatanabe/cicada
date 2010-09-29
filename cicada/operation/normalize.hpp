// -*- mode: c++ -*-

#ifndef __CICADA__OPERATION__NORMALIZE__HPP__
#define __CICADA__OPERATION__NORMALIZE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/generate.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

namespace cicada
{
  namespace operation
  {
    class Normalize : public Operation
    {
    public:
      Normalize(const std::string& parameter,
		const int __debug)
	: debug(__debug)
      {
	
      }
      

      template <typename FeaturePrefix, typename Feature>
      inline
      bool equal_prefix(const FeaturePrefix& prefix, const Feature& x)
      {
	return x.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), x.begin());
      }
      

      void operator()(data_type& data) const
      {
	typedef hypergraph_type::feature_set_type feature_set_type;

	hypergraph_type& hypergraph = data.hypergraph;

	if (debug)
	  std::cerr << "normalize" << std::endl;

	utils::resource start;
	
	hypergraph_type::node_set_type::const_iterator niter_end = hypergraph.nodes.end();
	for (hypergraph_type::node_set_type::const_iterator niter = hypergraph.nodes.begin(); niter != niter_end; ++ niter) {
	  const hypergraph_type::node_type& node = *niter;

	  hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	    hypergraph_type::edge_type& edge = hypergraph.edges[*eiter];

	    feature_set_type::data_type sum(0.0);
	    
	    feature_set_type::iterator fiter_end = edge.featuers.end();
	    for (feature_set_type::iterator fiter = edge.featuers.begin(); fiter != fiter_end; ++ fiter)
	      if (equal_prefix(feature_prefix, fiter->first))
		sum += fiter->second;
	    
	    if (sum != 0.0)
	      for (feature_set_type::iterator fiter = edge.featuers.begin(); fiter != fiter_end; ++ fiter)
		if (equal_prefix(feature_prefix, fiter->first))
		  fiter->second /= sum;
	  }
	}
	

	utils::resource end;
    
	if (debug)
	  std::cerr << "normalize cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << std::endl;
      }
    };
  };
};


#endif
