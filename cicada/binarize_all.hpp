// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_ALL__HPP__
#define __CICADA__BINARIZE_ALL__HPP__ 1

#include <cicada/binarize_base.hpp>

namespace cicada
{
  struct BinarizeAll : public BinarizeBase
  {
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, copy...
      target = source;
      
      phrase_type       binarized(2);
      hypergraph_type::edge_type::node_set_type tails(2);
      
      // we will traverse source-side in order to avoid confusion with newly created nodes...
      removed_type removed(source.edges.size(), false);

      position_set_type positions;
      node_chart_type   node_chart;
      label_chart_type  label_chart;
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  
	  if (edge_source.tails.size() <= 2) continue;
	  
	  removed[edge_source.id] = true;
	  
	  // we will create nodes in a chart structure, and exhaustively enumerate edges

	  const rule_type::symbol_set_type& rhs = edge_source.rule->rhs;
	  
	  // first, compute non-terminal spans...
	  positions.clear();
	  int pos = 0;
	  for (size_t i = 0; i != rhs.size(); ++ i)
	    if (rhs[i].is_non_terminal()) {
	      const symbol_type& non_terminal = rhs[i];
	      
	      // we will check non-terminal-index here...
	      const int non_terminal_index = non_terminal.non_terminal_index();
	      if (non_terminal_index > 0 && non_terminal_index - 1 != pos)
		throw std::runtime_error("hypergraph is not tail-sorted!");
	      ++ pos;
	      
	      positions.push_back(i);
	    }
	  
	  if (positions.size() != edge_source.tails.size())
	    throw std::runtime_error("invalid edge: # of non-terminals and tails size do not match");
	  
	  // seond, enumerate chart to compute node and edges...
	  node_chart.clear();
	  node_chart.resize(positions.size() + 1, hypergraph_type::invalid);

	  label_chart.clear();
	  label_chart.resize(positions.size() + 1);
	  
	  for (size_t i = 0; i != positions.size(); ++ i) {
	    node_chart(i, i + 1) = edge_source.tails[i];
	    label_chart(i, i + 1) = rhs[positions[i]].non_terminal();
	  }
	  
	  const symbol_type lhs_binarized = '[' + edge_source.rule->lhs.non_terminal_strip() + "^]";
	  
	  for (size_t length = 2; length <= positions.size(); ++ length)
	    for (size_t first = 0; first + length <= positions.size(); ++ first) {
	      const bool is_root = length == positions.size();
	      const size_t last = first + length;
	      
	      const hypergraph_type::id_type head = (is_root ? node_source.id : target.add_node().id);
	      const symbol_type& lhs = (is_root ? edge_source.rule->lhs : lhs_binarized);
	      
	      node_chart(first, last) = head;
	      label_chart(first, last) = lhs;
	      
	      for (size_t middle = first + 1; middle != last; ++ middle) {
		// [first, middle) and [middle, last)
		
		tails.front() = node_chart(first, middle);
		tails.back()  = node_chart(middle, last);
		
		binarized.clear();
		
		const size_t prefix_first = (is_root ? 0 : positions[first]);
		const size_t prefix_last  = positions[first];
		
		binarized.insert(binarized.end(), rhs.begin() + prefix_first, rhs.begin() + prefix_last);
		binarized.push_back(label_chart(first, middle));
		
		const size_t middle_first = positions[middle - 1] + 1;
		const size_t middle_last  = positions[middle];
		
		binarized.insert(binarized.end(), rhs.begin() + middle_first, rhs.begin() + middle_last);
		binarized.push_back(label_chart(middle, last));
		
		const size_t suffix_first = positions[last - 1] + 1;
		const size_t suffix_last  = (is_root ? static_cast<int>(rhs.size()) : positions[last - 1] + 1);
		
		binarized.insert(binarized.end(), rhs.begin() + suffix_first, rhs.begin() + suffix_last);
		
		hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
		edge_new.rule = rule_type::create(rule_type(lhs, binarized.begin(), binarized.end()));
		
		target.connect_edge(edge_new.id, head);

		// assign features...
		if (is_root) {
		  edge_new.features   = edge_source.features;
		  edge_new.attributes = edge_source.attributes;
		}
	      }
	    }
	}
      }
      
      // further resize...
      removed.resize(target.edges.size(), false);
      
      hypergraph_type graph_removed;
      
      topologically_sort(target, graph_removed, filter(removed));
      
      target.swap(graph_removed);
    }
  };

  inline
  void binarize_all(const HyperGraph& source, HyperGraph& target)
  {
    BinarizeAll binarizer;
    
    binarizer(source, target);
  }

  inline
  void binarize_all(HyperGraph& source)
  {
    HyperGraph target;

    BinarizeAll binarizer;
    
    binarizer(source, target);
    
    source.swap(target);
  }  
};

#endif
