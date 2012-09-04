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
    
    typedef std::vector<id_type, std::allocator<id_type> > edge_list_type;
    typedef std::vector<edge_list_type, std::allocator<edge_list_type> > edge_set_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;

      edge_set_type outgoing(target.nodes.size());
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = target.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = target.edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = *eiter;
	
	edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	for (edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	  outgoing[*titer].push_back(edge.id);
      }
      
      hypergraph_type::node_set_type::const_reverse_iterator niter_end = target.nodes.rend();
      for (hypergraph_type::node_set_type::const_reverse_iterator niter = target.nodes.rbegin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	  
	// this should not happen, though..
	if (node.edges.empty()) continue;
	if (outgoing[node.id].empty()) continue;
	  
	if (outgoing[node.id].size() == 1) {

	  if (! target.edges[outgoing[node.id].front()].features.empty()) {
	    const feature_set_type intersected(target.edges[outgoing[node.id].front()].features);
	      
	    node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	      target.edges[*eiter].features += intersected;
	      
	    target.edges[outgoing[node.id].front()].features.clear();
	  }
	} else {
	  feature_set_type intersected(target.edges[outgoing[node.id].front()].features);
	    
	  edge_list_type::const_iterator oiter_end = outgoing[node.id].end();
	  for (edge_list_type::const_iterator oiter = outgoing[node.id].begin() + 1; oiter != oiter_end; ++ oiter)
	    intersected.intersect(target.edges[*oiter].features);

	  if (! intersected.empty()) {
	    edge_list_type::const_iterator oiter_end = outgoing[node.id].end();
	    for (edge_list_type::const_iterator oiter = outgoing[node.id].begin(); oiter != oiter_end; ++ oiter)
	      target.edges[*oiter].features -= intersected;
	      
	    node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	      target.edges[*eiter].features += intersected;
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
