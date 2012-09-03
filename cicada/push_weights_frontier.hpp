// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PUSH_WEIGHTS_FRONTIER__HPP__
#define __CICADA__PUSH_WEIGHTS_FRONTIER__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>

namespace cicada
{
  
  struct PushWeightsFrontier
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    typedef std::vector<feature_set_type, std::allocator<feature_set_type> > feature_map_type;
    
    typedef std::vector<bool, std::allocator<bool> > processed_type;

    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;

      feature_map_type outside(target.nodes.size());
      processed_type   processed(target.nodes.size());
      
      // visit in a topological order, then, compute intersection of features, 
      processed[target.goal] = true;
      
      {
	hypergraph_type::node_set_type::const_reverse_iterator niter_end = target.nodes.rend();
	for (hypergraph_type::node_set_type::const_reverse_iterator niter = target.nodes.rbegin(); niter != niter_end; ++ niter) {
	  const hypergraph_type::node_type& node = *niter;
	  
	  // this should not happen, though..
	  if (node.edges.empty()) continue;
	  
	  node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	    edge_type& edge = target.edges[*eiter];
	    
	    feature_set_type intersected(edge.features);
	    intersected += outside[node.id];
	    
	    edge_type::node_set_type::const_iterator titer_begin = edge.tails.begin();
	    edge_type::node_set_type::const_iterator titer_end   = edge.tails.end();
	    for (edge_type::node_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	      if (! processed[*titer]) {
		outside[*titer] = intersected;
		processed[*titer] = true;
	      } else
		outside[*titer].intersect(intersected);
	    }
	  }
	}
      }
    }
  };
  
  inline
  void push_weights_frontier(const HyperGraph& source, HyperGraph& target)
  {
    PushWeightsFrontier()(source, target);
  }
  
  inline
  void push_weights_frontier(HyperGraph& graph)
  {
    HyperGraph x;
    push_weights_frontier(graph, x);
    graph.swap(x);
  }
};

#endif
