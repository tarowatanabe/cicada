// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SORT_TAIL__HPP__
#define __CICADA__SORT_TAIL__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>

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
      typedef std::vector<id_type, std::allocator<id_type> >         tails_type;
      typedef std::vector<int, std::allocator<int> >                 index_type;

      graph = x;
      
      if (! graph.is_valid()) return;
      
      rhs_type   rhs;
      tails_type tails;
      index_type index;
      
      hypergraph_type::edge_set_type::iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	edge_type& edge = *eiter;

	if (edge.tails.size() <= 1) continue;
	
	rhs.clear();
	tails.clear();
	index.clear();
	
	const symbol_type& lhs = edge.rule->lhs;
	rhs.insert(rhs.end(), edge.rule->rhs.begin(), edge.rule->rhs.end());
	tails.insert(tails.end(), edge.tails.begin(), edge.tails.end());
	index.resize(tails.size());
	
	int pos = 0;
	rhs_type::iterator riter_end = rhs.end();
	for (rhs_type::iterator riter = rhs.begin(); riter != riter_end; ++ riter)
	  if (riter->is_non_terminal()) {
	    const int non_terminal_pos = riter->non_terminal_index();
	    
	    index[non_terminal_pos == 0 ? pos : non_terminal_pos - 1] = pos;
	    
	    *riter = riter->non_terminal(pos + 1);
	    ++ pos;
	  }
	
	// sort tail wrt index...
	for (size_t i = 0; i != tails.size(); ++ i)
	  tails[i] = edge.tails[index[i]];
	
	// re-assign edge's rule and edge's tails...
	
	edge.rule = rule_type::create(rule_type(lhs, rhs.begin(), rhs.end()));
	edge.tails.assign(tails.begin(), tails.end());
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
