// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SPAN_FOREST__HPP__
#define __CICADA__SPAN_FOREST__HPP__ 1

#include <stdint.h>

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/vocab.hpp>
#include <cicada/symbol.hpp>

namespace cicada
{
  
  struct SpanForest
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    typedef hypergraph_type::attribute_set_type attribute_set_type;
    
    typedef attribute_set_type::attribute_type attribute_type;

    typedef Vocab vocab_type;

    typedef std::pair<int, int> span_type;
    typedef std::vector<span_type, std::allocator<span_type> > span_set_type;

    typedef std::vector<bool, std::allocator<bool> > visited_type;

    SpanForest() : attr_span_first("span-first"), attr_span_last("span-last") {}
    
    void operator()(const hypergraph_type& __graph, hypergraph_type& graph)
    {
      graph = __graph;
      if (graph.is_valid()) {
	visited.clear();
	visited.resize(graph.nodes.size(), false);
	
	spans_node.clear();
	spans_node.reserve(graph.nodes.size());
	spans_node.resize(graph.nodes.size(), span_type(0, 0));
	
	traverse(graph, graph.goal);
      }
    }
    
    // top-down traversal to compute spans...
    void traverse(hypergraph_type& graph, id_type node_id)
    {
      if (visited[node_id]) return;

      const node_type& node = graph.nodes[node_id];
      
      // parant-span starts from spans_node[node_id].first;
      
      node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	edge_type& edge = graph.edges[*eiter];
	
	int span_pos = spans_node[node_id].first;
	
	edge.attributes[attr_span_first] = attribute_set_type::int_type(span_pos);
	
	int non_terminal_pos = 0;
	rule_type::symbol_set_type::const_iterator siter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator siter = edge.rule->rhs.begin(); siter != siter_end; ++ siter) {
	  if (siter->is_non_terminal()) {
	    const int __non_terminal_index = siter->non_terminal_index();
	    const int non_terminal_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	    
	    spans_node[edge.tails[non_terminal_index]].first = span_pos;
	    
	    traverse(graph, edge.tails[non_terminal_index]);
	    
	    span_pos = spans_node[edge.tails[non_terminal_index]].second;
	    
	    ++ non_terminal_pos;
	  } else if (*siter != vocab_type::EPSILON)
	    ++ span_pos;
	}
	
	edge.attributes[attr_span_last] = attribute_set_type::int_type(span_pos);
	
	spans_node[node_id].second = utils::bithack::max(spans_node[node_id].second, span_pos);
      }

      visited[node_id] = true;
    }

    span_set_type spans_node;
    visited_type  visited;
    
    const attribute_type attr_span_first;
    const attribute_type attr_span_last;
  };

  inline
  void span_forest(HyperGraph& graph)
  {
    HyperGraph target;
    SpanForest __span_forest;
    
    __span_forest(graph, target);
    
    graph.swap(target);
  }
  
  inline
  void span_forest(const HyperGraph& graph, HyperGraph& target)
  {
    SpanForest __span_forest;
    
    __span_forest(graph, target);
  }
};

#endif
