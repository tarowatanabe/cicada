// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_CYK__HPP__
#define __CICADA__BINARIZE_CYK__HPP__ 1

#include <deque>

#include <cicada/binarize_base.hpp>

namespace cicada
{
  struct BinarizeCYK : public BinarizeBase
  {
    BinarizeCYK(const int __order=0)
      : order(__order) {}

    typedef std::vector<bool, std::allocator<bool> > visited_type;
    
    typedef std::pair<int, int> span_type;
    typedef std::vector<span_type, std::allocator<span_type> > span_set_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > label_set_type;
    typedef std::deque<label_set_type, std::allocator<label_set_type> > label_nodes_type;
    
    typedef utils::chart<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_chart_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, topological-order traversal to compute span for each node and 
      // ancestors.
      // I don't know how to handle unary rules... we will take label from the bottom. (top or concatenate them together?)
      
      // order <= 0 implies all the ancestors...

      // initialize target-nodes...
      target.clear();
      
      if (! source.is_valid()) return;
      
      target.nodes.resize(source.nodes.size());
      for (size_t i = 0; i != target.nodes.size(); ++ i)
	target.nodes[i].id = i;
      target.goal = source.goal;
      
      visited.clear();
      visited.reserve(source.nodes.size());
      visited.resize(source.nodes.size(), false);
      
      spans.clear();
      spans.reserve(source.nodes.size());
      spans.resize(source.nodes.size(), span_type(0, 0));
      
      // annotate spans, label etc.
      traverse_span(source, source.goal);
      
    }

    void traverse_span(const hypergraph_type& graph, const hypergraph_type::id_type& node_id)
    {
      if (visited[node_id]) return;
      
      const hypergraph_type::node_type& node = graph.nodes[node_id];
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	
	int span_pos = spans[node_id].first;
	
	int non_terminal_pos = 0;
	rule_type::symbol_set_type::const_iterator siter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator siter = edge.rule->rhs.begin(); siter != siter_end; ++ siter) {
	  if (siter->is_non_terminal()) {
	    int non_terminal_index = siter->non_terminal_index() - 1;
	    if (non_terminal_index < 0)
	      non_terminal_index = non_terminal_pos;
	    ++ non_terminal_pos;
	    
	    spans[edge.tails[non_terminal_index]].first = span_pos;
	    
	    traverse(graph, edge.tails[non_terminal_index]);
	    
	    span_pos = spans[edge.tails[non_terminal_index]].second;
	  } else if (*siter != vocab_type::EPSILON)
	    ++ span_pos;
	}
	
	spans[node_id].second = utils::bithack::max(spans[node_id].second, span_pos);
      }
      
      visited[node_id] = true;
    }
	

    visited_type  visited;
    span_set_type spans;
    
    const int order;
  };
  
  inline
  void binarize_cyk(const HyperGraph& source, HyperGraph& target, const int order=0)
  {
    BinarizeCYK binarizer(order);
    
    binarizer(source, target);
  }

  inline
  void binarize_cyk(HyperGraph& source, const int order=0)
  {
    HyperGraph target;

    BinarizeCYK binarizer(order);
    
    binarizer(source, target);
    
    source.swap(target);
  }
};

#endif
