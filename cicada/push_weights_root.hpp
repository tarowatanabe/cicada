// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PUSH_WEIGHTS_ROOT__HPP__
#define __CICADA__PUSH_WEIGHTS_ROOT__HPP__ 1

#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <utility>

#include <cicada/hypergraph.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sort_topologically.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/mathop.hpp>
#include <utils/bithack.hpp>


namespace cicada
{
  
  struct PushWeightsRoot
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    typedef std::vector<feature_set_type, std::allocator<feature_set_type> > feature_map_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;

      feature_map_type features(target.nodes.size());
      
      // visit in a topological order, then, computed intersection of features, 
      hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	// this should not happen, though..
	if (node.edges.empty()) continue;
	
	// collect pushed features...
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  edge_type& edge = target.edges[*eiter];
	  
	  edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	  for (edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	    edge.features += features[*titer];
	}
	
	if (node.id == target.goal) continue;
	
	// push features...
	
	if (node.edges.size() == 1) {
	  features[node.id] += target.edges[node.edges.front()].features;
	  
	  target.edges[node.edges.front()].features.clear();
	} else {
	  feature_set_type intersected(target.edges[node.edges.front()].features);
	  
	  node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (node_type::edge_set_type::const_iterator eiter = node.edges.begin() + 1; eiter != eiter_end; ++ eiter)
	    intersected.intersect(target.edges[*eiter].features);
	  
	  if (! intersected.empty()) {
	    features[node.id] += intersected;
	    
	    node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	      target.edges[*eiter].features -= intersected;
	  }
	}
      }
    }
  };
  
  
  inline
  void push_weights_root(const HyperGraph& source, HyperGraph& target)
  {
    PushWeightsRoot()(source, target);
  }
  
  inline
  void push_weights_root(HyperGraph& graph)
  {
    HyperGraph x;
    push_weights_root(graph, x);
    graph.swap(x);
  }
};

#endif
