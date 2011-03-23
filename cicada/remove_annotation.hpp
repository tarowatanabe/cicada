// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__REMOVE_ANNOTATION__HPP__
#define __CICADA__REMOVE_ANNOTATION__HPP__ 1

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/vocab.hpp>
#include <cicada/symbol.hpp>

namespace cicada
{
  struct RemoveAnnotation
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    typedef rule_type::symbol_type     symbol_type;
    typedef rule_type::symbol_set_type symbol_set_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      target = source;
      
      hypergraph_type::edge_set_type::iterator eiter_end = target.edges.end();
      for (hypergraph_type::edge_set_type::iterator eiter = target.edges.begin(); eiter != eiter_end; ++ eiter) {
	edge_type& edge = *eiter;
	const rule_type& rule = *edge.rule;
	
	const symbol_type lhs = remove_annotation(symbol);
	symbol_set_type rhs(rule.rhs);
	
	symbol_set_type::iterator riter_end = rhs.end();
	for (symbol_set_type::iterator riter = rhs.begin(); riter != riter_end; ++ riter)
	  if (riter->is_non_terminal())
	    *riter = remove_annotation(*riter);
	
	edge.rule = rule_type::create(rule_type(lhs, rhs));
      }
    }
    
    symbol_type remove_annotation(const symbol_type& symbol)
    {
      const utils::piece piece = symbol.non_terminal().non_terminal_strip();
      
      utils::piece::size_type pos = piece.find('@');
      if (pos != utils::piece::npos())
	return '[' + piece.substr(0, pos) + ']';
      else
	return symbol.non_terminal();
    }
    
  };
};

#endif
