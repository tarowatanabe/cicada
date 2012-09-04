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
      edge_set_type leaning(target.nodes.size());
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = target.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = target.edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = *eiter;
	
	if (edge.tails.empty()) continue;
	
	int tail_pos = 0;
	rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) 
	  if (riter->is_non_terminal()) {
	    const int __non_terminal_index = riter->non_terminal_index();
	    const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, tail_pos, __non_terminal_index - 1);
	    
	    if (antecedent_index == 0)
	      leaning[edge.tails[tail_pos]].push_back(edge.id);
	    
	    outgoing[edge.tails[tail_pos]].push_back(edge.id);
	    
	    ++ tail_pos;
	  }
      }
      
      hypergraph_type::node_set_type::const_reverse_iterator niter_end = target.nodes.rend();
      for (hypergraph_type::node_set_type::const_reverse_iterator niter = target.nodes.rbegin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	// this should not happen, though..
	if (node.edges.empty()) continue;
	if (outgoing[node.id].empty()) continue;
	
	const edge_list_type& accumulate = (leaning[node.id].empty() ? outgoing[node.id] : leaning[node.id]);
	
	if (accumulate.size() == 1) {
	  
	  if (! target.edges[accumulate.front()].features.empty()) {
	    const feature_set_type intersected(target.edges[accumulate.front()].features);
	    
	    node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	      target.edges[*eiter].features += intersected;
	    
	    edge_list_type::const_iterator oiter_end = outgoing[node.id].end();
	    for (edge_list_type::const_iterator oiter = outgoing[node.id].begin(); oiter != oiter_end; ++ oiter)
	      target.edges[*oiter].features -= intersected;
	  }
	} else {
	  feature_set_type intersected(target.edges[accumulate.front()].features);
	    
	  edge_list_type::const_iterator oiter_end = accumulate.end();
	  for (edge_list_type::const_iterator oiter = accumulate.begin() + 1; oiter != oiter_end; ++ oiter)
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
