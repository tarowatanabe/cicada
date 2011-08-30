// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SORT_TAIL__HPP__
#define __CICADA__SORT_TAIL__HPP__ 1

#include <vector>
#include <algorithm>

#include <cicada/hypergraph.hpp>

#include <utils/bithack.hpp>

namespace cicada
{
  struct SortTail
  {
    typedef HyperGraph hypergraph_type;

    
    typedef hypergraph_type::id_type     id_type;
    typedef hypergraph_type::symbol_type symbol_type;
    typedef hypergraph_type::rule_type   rule_type;
    typedef hypergraph_type::edge_type   edge_type;
    typedef hypergraph_type::node_type   node_type;
    
    void operator()(const hypergraph_type& x, hypergraph_type& graph)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > rhs_type;

      graph = x;
      
      if (! graph.is_valid()) return;
      
      rhs_type   rhs;
      
      hypergraph_type::edge_set_type::iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	edge_type& edge = *eiter;

	if (edge.tails.empty()) continue;
	
	const symbol_type& lhs = edge.rule->lhs;

	rhs.clear();
	rhs.insert(rhs.end(), edge.rule->rhs.begin(), edge.rule->rhs.end());
	
	edge_type::node_set_type tails(edge.tails.size());
	
	int pos = 0;
	rhs_type::iterator riter_end = rhs.end();
	for (rhs_type::iterator riter = rhs.begin(); riter != riter_end; ++ riter)
	  if (riter->is_non_terminal()) {
	    const int non_terminal_pos = riter->non_terminal_index();
	    const int antecedent_index = utils::bithack::branch(non_terminal_pos == 0, pos, non_terminal_pos - 1);
	    
	    tails[pos] = edge.tails[antecedent_index];
	    
	    *riter = riter->non_terminal();
	    ++ pos;
	  }
	
	edge.rule = rule_type::create(rule_type(lhs, rhs.begin(), rhs.end()));
	edge.tails = tails;
      }
    }    
  };

  inline
  void sort_tail(const HyperGraph& x, HyperGraph& y)
  {
    SortTail()(x, y);
  }
    
  inline
  void sort_tail(HyperGraph& graph)
  {
    HyperGraph x;
    sort_tail(graph, x);
    graph.swap(x);
  }

};

#endif
