// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TREE_GRAMMAR_SIMPLE__HPP__
#define __CICADA__TREE_GRAMMAR_SIMPLE__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/tree_grammar_mutable.hpp>

#include <google/dense_hash_set>

namespace cicada
{
  class TreeGrammarFallback : public TreeGrammarMutable
  {
  public:
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::rule_type     graph_rule_type;
    typedef hypergraph_type::rule_ptr_type graph_rule_ptr_type;

  public:
    struct rule_ptr_hash
    {
      size_t operator()(const graph_rule_ptr_type& x) const
      {
	return (x ? hash_value(*x) : size_t(0));
      }
    };

    struct rule_ptr_equal
    {
      bool operator()(const graph_rule_ptr_type& x, const graph_rule_ptr_type& y) const
      {
	return (x == y || x && y && *x == *y);
      }
    };

    typedef google::dense_hash_set<graph_rule_ptr_type, rule_ptr_hash, rule_ptr_equal> graph_rule_ptr_set_type;

  public:

    // copy, but use default non-terminal category
    TreeGrammarFallback(const hypergraph_type& graph, const symbol_type& non_terminal)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;

      graph_rule_ptr_set_type rules;
      rules.set_empty_key(graph_rule_ptr_type());

      feature_set_type features;
      features["tree-insertion-penalty"] = -1.0;

      non_terminal_set_type non_terminals;
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = *eiter;
	
	if (rules.find(edge.rule) != rules.end()) continue;
	rules.insert(edge.rule);
	
	non_terminals.clear();
	symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	for (symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
	  non_terminals.push_back(riter->is_non_terminal() ? non_terminal.non_terminal(riter->non_terminal_index()) : *riter);
	
	rule_ptr_type rule_source(rule_type::create(rule_type(edge.rule->lhs, edge.rule->rhs.begin(), edge.rule->rhs.end())));
	rule_ptr_type rule_target(rule_type::create(rule_type(non_terminal, non_terminals.begin(), non_terminals.end())));
	
	insert(rule_pair_type(rule_source, rule_target, features));
      }
    }
    
    // simply copy and preserve the same non-terminal category...
    TreeGrammarFallback(const hypergraph_type& graph)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;

      graph_rule_ptr_set_type rules;
      rules.set_empty_key(graph_rule_ptr_type());

      feature_set_type features;
      features["tree-insertion-penalty"] = -1.0;

      non_terminal_set_type non_terminals;
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = *eiter;
	
	if (rules.find(edge.rule) != rules.end()) continue;
	rules.insert(edge.rule);
	
	rule_ptr_type rule(rule_type::create(rule_type(edge.rule->lhs, edge.rule->rhs.begin(), edge.rule->rhs.end())));
	
	insert(rule_pair_type(rule, rule, features));
      }
    }

  };
  
};

#endif
