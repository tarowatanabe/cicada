// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PUSH_WEIGHTS_LEFT__HPP__
#define __CICADA__PUSH_WEIGHTS_LEFT__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>

#include <utils/bithack.hpp>

namespace cicada
{
  
  struct PushWeightsLeft
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef std::vector<id_type, std::allocator<id_type> > edge_set_type;
    typedef std::vector<edge_set_type, std::allocator<edge_set_type> > edge_map_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;
      
      edge_map_type    edges(source.nodes.size());

      // visit in a topological order, then, compute a "left-learning graph"
      {
	hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
	for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	  const hypergraph_type::node_type& node = *niter;
	
	  node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	    const edge_type& edge = target.edges[*eiter];
	  
	    if (edge.tails.empty()) continue;
	  
	    int tail_pos = 0;
	    rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	    for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) 
	      if (riter->is_non_terminal()) {
		const int __non_terminal_index = riter->non_terminal_index();
		const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, tail_pos, __non_terminal_index - 1);
		
		if (antecedent_index == 0) break;
		
		++ tail_pos;
	      }
	    
	    edges[edge.tails[tail_pos]].push_back(edge.id);
	  }
	}
      }
      
      // visit in an inverse topological order, and push-features!
      hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	// nothing to propagate!
	if (edges[node.id].empty()) continue;
	
	if (edges[node.id].size() == 1) {
	  node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	    target.edges[*eiter].features += target.edges[edges[node.id].front()].features;
	  
	  target.edges[edges[node.id].front()].features.clear();
	} else {
	  feature_set_type intersected(target.edges[edges[node.id].front()].features);
	  
	  edge_set_type::const_iterator eiter_end = edges[node.id].end();
	  for (edge_set_type::const_iterator eiter = edges[node.id].begin() + 1; eiter != eiter_end; ++ eiter)
	    intersected.intersect(target.edges[*eiter].features);
	  
	  if (! intersected.empty()) {
	    edge_set_type::const_iterator eiter_end = edges[node.id].end();
	    for (edge_set_type::const_iterator eiter = edges[node.id].begin(); eiter != eiter_end; ++ eiter)
	      target.edges[*eiter].features -= intersected;
	    
	    {
	      node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	      for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
		target.edges[*eiter].features += intersected;
	    }
	  }
	}
      }
    }
  };
  
  
  inline
  void push_weights_left(const HyperGraph& source, HyperGraph& target)
  {
    PushWeightsLeft()(source, target);
  }
  
  inline
  void push_weights_left(HyperGraph& graph)
  {
    HyperGraph x;
    push_weights_left(graph, x);
    graph.swap(x);
  }
};

#endif
