// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PERMUTE_DETERMINISTIC__HPP__
#define __CICADA__PERMUTE_DETERMINISTIC__HPP__ 1

#include <algorithm>

#include <cicada/hypergraph.hpp>

#include <utils/bithack.hpp>


namespace cicada
{
  // deterministically permute
  
  template <typename Filter>
  struct PermuteDeterministic
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::symbol_type      symbol_type;
    
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;

    typedef Filter filter_type;
    
    PermuteDeterministic(const filter_type& __filter)
      : filter(__filter) {}
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      target = source;
      
      if (! target.is_valid()) return;
      
      hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  hypergraph_type::edge_type& edge = target.edges[*eiter];
	  
	  if (! filter(edge.rule->lhs)) continue;
	  
	  rule_type::symbol_set_type rhs(edge.rule->rhs);
	  
	  // assign index...
	  int non_terminal_pos = 0;
	  rule_type::symbol_set_type::iterator riter_end = rhs.end();
	  for (rule_type::symbol_set_type::iterator riter = rhs.begin(); riter != riter_end; ++ riter)
	    if (riter->is_non_terminal()) {
	      const int __non_terminal_index = riter->non_terminal_index();
	      const int pos = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	      ++ non_terminal_pos;
	      
	      *riter = riter->non_terminal(pos + 1);
	    }
	  
	  // inverse
	  std::reverse(rhs.begin(), rhs.end());
	  
	  // new rule!
	  edge.rule = rule_type::create(rule_type(edge.rule->lhs, rhs));
	}
      }
    }
    
  private:
    const filter_type& filter;
  };
  
  template <typename Filter>
  inline
  void permute_deterministic(const HyperGraph& source, HyperGraph& target, const Filter& filter)
  {
    PermuteDeterministic<Filter> permutation(filter);
    
    permutation(source, target);
  }
  
  template <typename Filter>
  inline
  void permute_deterministic(HyperGraph& graph, const Filter& filter)
  {
    PermuteDeterministic<Filter> permutation(filter);
    
    HyperGraph permuted;
    permutation(graph, permuted);
    permuted.swap(graph);
  }

};


#endif
