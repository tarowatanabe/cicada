// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_RIGHT__HPP__
#define __CICADA__BINARIZE_RIGHT__HPP__ 1

#include <cicada/binarize_base.hpp>

namespace cicada
{
  struct BinarizeRight : public BinarizeBase
  {
    BinarizeRight(const int __order=-1)
      : order(__order) {}

    template <typename Iterator>
    symbol_type binarized_label(const symbol_type& lhs, Iterator first, Iterator last)
    {
      std::string label = lhs.non_terminal().non_terminal_strip() + '^';
      for (/**/; first != last; ++ first) {
	label += '_';
	label += first->non_terminal().non_terminal_strip();
      }
      return '[' + label + ']';
    }
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, copy...
      target = source;
      
      context_type context;
      phrase_type binarized(2);
      hypergraph_type::edge_type::node_set_type tails(2);
      
      // we will traverse source-side in order to avoid confusion with newly created nodes...
      removed_type removed(source.edges.size(), false);
      
      position_set_type positions;
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  	  
	  if (edge_source.tails.size() <= 2) continue;
	  
	  removed[edge_source.id] = true;
	  
	  // right most antecedents binarization...
	  
	  const symbol_type& lhs = edge_source.rule->lhs;
	  const rule_type::symbol_set_type& rhs = edge_source.rule->rhs;
	  
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
	  
	  hypergraph_type::id_type head = edge_source.head;
	  symbol_type              head_label = lhs;
	  
	  
	  context.clear();
	  
	  for (size_t length = positions.size(); length >= 2; -- length) {
	    const bool is_root = length == positions.size();
	    const size_t last = positions.size();
	    const size_t first = last - length;
	    const size_t middle = first + 1;
	    
	    context.push_back(rhs[positions[first]].non_terminal());
	    
	    const symbol_type non_terminal = (length == 2
					      ? rhs[positions[last - 1]].non_terminal()
					      : (order < 0
						 ? binarized_label(lhs, context.begin(), context.end())
						 : binarized_label(lhs, std::max(context.begin(), context.end() - order), context.end())));
	    
	    binarized.clear();
	    
	    const size_t prefix_first = (is_root ? 0 : positions[first]);
	    const size_t prefix_last  = positions[first];
	    
	    binarized.insert(binarized.end(), rhs.begin() + prefix_first, rhs.begin() + prefix_last);
	    binarized.push_back(rhs[positions[first]].non_terminal());
	    
	    const size_t middle_first = positions[middle - 1] + 1;
	    const size_t middle_last  = positions[middle];
	    
	    binarized.insert(binarized.end(), rhs.begin() + middle_first, rhs.begin() + middle_last);
	    binarized.push_back(non_terminal);
	    
	    const size_t suffix_first = positions[last - 1] + 1;
	    const size_t suffix_last  = (is_root ? static_cast<int>(rhs.size()) : positions[last - 1] + 1);
	    
	    binarized.insert(binarized.end(), rhs.begin() + suffix_first, rhs.begin() + suffix_last);
	    
	    tails.front() = edge_source.tails[first];
	    tails.back() = length == 2 ? edge_source.tails[last - 1] : target.add_node().id;
	    
	    hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	    
	    edge_new.rule = rule_type::create(rule_type(head_label, binarized.begin(), binarized.end()));
	    
	    if (is_root) {
	      edge_new.features   = edge_source.features;
	      edge_new.attributes = edge_source.attributes;
	    }
	    
	    target.connect_edge(edge_new.id, head);
	    
	    head = tails.back();
	    head_label = non_terminal;
	  }
	}
      }
      
      // further resize...
      removed.resize(target.edges.size(), false);
      
      hypergraph_type graph_removed;
      
      topologically_sort(target, graph_removed, filter(removed));
      
      target.swap(graph_removed);
    }
    
    int order;
  };
  
  inline
  void binarize_right(const HyperGraph& source, HyperGraph& target, const int order=-1)
  {
    BinarizeRight binarizer(order);
    
    binarizer(source, target);
  }

  inline
  void binarize_right(HyperGraph& source, const int order=-1)
  {
    HyperGraph target;

    BinarizeRight binarizer(order);
    
    binarizer(source, target);
    
    source.swap(target);
  }

};

#endif
