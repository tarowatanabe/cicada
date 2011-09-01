// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__POSTERIOR__HPP__
#define __CICADA__POSTERIOR__HPP__ 1

#include <cicada/semiring.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/inside_outside.hpp>

namespace cicada
{
  template <typename Semiring, typename Function>
  struct Posterior
  {
    // compute posterior feature...
    // we will run inside/outside and re-assign features for each hyperedge by clearning all the
    // features and assign "posterior" feature

  public:
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::feature_set_type feature_set_type;
    typedef feature_set_type::feature_type    feature_type;
    
    typedef Function function_type;
    typedef Semiring weight_type;
    
    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
    
    Posterior(Function __function) : function(__function), feat_posterior("posterior") {}
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      target = source;
      
      if (! target.is_valid()) return;

      hypergraph_type& graph = target;
      
      inside.clear();
      inside.reserve(graph.nodes.size());
      inside.resize(graph.nodes.size());
      
      outside.clear();
      outside.reserve(graph.nodes.size());
      outside.resize(graph.nodes.size());
      
      cicada::inside(graph, inside, function);
      cicada::outside(graph, inside, outside, function);
      
      // reassign features...
      const weight_type weight_total = inside.back();

      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {

	const weight_type weight_outside = outside[niter->id];
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = niter->edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = niter->edges.begin(); eiter != eiter_end; ++ eiter) {
	  hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  weight_type weight = weight_outside * function(edge) / weight_total;
	  hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	  for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	    weight *= inside[*titer];
	  
	  // assign features!
	  edge.features.clear();
	  
	  const double logweight = cicada::semiring::log(weight);
	  
	  if (logweight != 0.0)
	    edge.features[feat_posterior] = logweight;
	}
      }
    }
    
    function_type function;
    
    weight_set_type inside;
    weight_set_type outside;
    
    const feature_type feat_posterior;
  };
  
  
  template <typename Function>
  inline
  void posterior(const HyperGraph& source, HyperGraph& target, const Function& func)
  {
    Posterior<typename Function::value_type, Function> __posterior(func);
    __posterior(source, target);
  }
  
  template <typename Function>
  inline
  void posterior(HyperGraph& graph, const Function& func)
  {
    Posterior<typename Function::value_type, Function> __posterior(func);
    HyperGraph computed;
    __posterior(graph, computed);
    graph.swap(computed);
  }

};

#endif
