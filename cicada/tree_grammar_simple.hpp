// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TREE_GRAMMAR_SIMPLE__HPP__
#define __CICADA__TREE_GRAMMAR_SIMPLE__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/tree_grammar_mutable.hpp>

#include <utils/unordered_set.hpp>

namespace cicada
{
  class TreeGrammarGlue : public TreeGrammarMutable
  {
    
  public:
    TreeGrammarGlue(const symbol_type& goal_source,
		    const symbol_type& goal_target,
		    const symbol_type& non_terminal_source,
		    const symbol_type& non_terminal_target,
		    const bool __straight,
		    const bool __inverted)
      : TreeGrammarMutable(true),
	straight(__straight),
	inverted(__inverted)
    {
      construct(goal_source, goal_target, non_terminal_source, non_terminal_target);
    }
    
    transducer_ptr_type clone() const { return transducer_ptr_type(new TreeGrammarGlue(*this)); }

    bool valid_span(int first, int last, int distance) const
    {
      return (straight && inverted ? true : first == 0);
    }

  private:
    void construct(const symbol_type& goal_source,
		   const symbol_type& goal_target,
		   const symbol_type& non_terminal_source,
		   const symbol_type& non_terminal_target)
    {
      if (! goal_source.is_non_terminal())
	throw std::runtime_error("invalid non-terminal: " + static_cast<const std::string&>(goal_source));
      if (! goal_target.is_non_terminal())
	throw std::runtime_error("invalid non-terminal: " + static_cast<const std::string&>(goal_target));
      
      if (! non_terminal_source.is_non_terminal())
	throw std::runtime_error("invalid non-terminal: " + static_cast<const std::string&>(non_terminal_source));
      if (! non_terminal_target.is_non_terminal())
	throw std::runtime_error("invalid non-terminal: " + static_cast<const std::string&>(non_terminal_target));
      
      feature_set_type features;
      features["tree-glue-unary-penalty"] = -1;

      attribute_set_type attributes;
      attributes["tree-fallback"] = attribute_set_type::int_type(1);
      
      // unary rule...
      insert(rule_pair_type(rule_type::create(rule_type(goal_source, &non_terminal_source, (&non_terminal_source) + 1)),
			    rule_type::create(rule_type(goal_target, &non_terminal_target, (&non_terminal_target) + 1)),
			    features,
			    attributes));

      if (straight) {
	std::vector<symbol_type, std::allocator<symbol_type> > source(2);
	std::vector<symbol_type, std::allocator<symbol_type> > target(2);
	
	source.front() = goal_source.non_terminal(1);
	source.back()  = non_terminal_source.non_terminal(2);

	target.front() = goal_target.non_terminal(1);
	target.back()  = non_terminal_target.non_terminal(2);
	
	feature_set_type features;
	features["tree-glue-straight-penalty"] = -1;

	insert(rule_pair_type(rule_type::create(rule_type(goal_source, source.begin(), source.end())),
			      rule_type::create(rule_type(goal_target, target.begin(), target.end())),
			      features,
			      attributes));
      }
      
      if (inverted) {
	std::vector<symbol_type, std::allocator<symbol_type> > source(2);
	std::vector<symbol_type, std::allocator<symbol_type> > target(2);
	
	source.front() = goal_source.non_terminal(1);
	source.back()  = non_terminal_source.non_terminal(2);
	
	target.front() = non_terminal_target.non_terminal(2);
	target.back()  = goal_target.non_terminal(1);
	
	feature_set_type features;
	features["tree-glue-inverted-penalty"] = -1;

	insert(rule_pair_type(rule_type::create(rule_type(goal_source, source.begin(), source.end())),
			      rule_type::create(rule_type(goal_target, target.begin(), target.end())),
			      features,
			      attributes));	
      }
    }
    
  private:
    bool straight;
    bool inverted;
  };

  class TreeGrammarFallback : public TreeGrammarMutable
  {
  public:
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::rule_type     graph_rule_type;
    typedef hypergraph_type::rule_ptr_type graph_rule_ptr_type;

    typedef feature_set_type::feature_type feature_type;

  private:
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
	return (x == y || (x && y && *x == *y));
      }
    };

    typedef utils::unordered_set<graph_rule_ptr_type, rule_ptr_hash, rule_ptr_equal, std::allocator<graph_rule_ptr_type> >::type graph_rule_ptr_set_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
  public:

    TreeGrammarFallback(const symbol_type& __goal, const symbol_type& __non_terminal)
      : goal(__goal), non_terminal(__non_terminal),
	feat_penalty("tree-insertion-penalty"),
	feat_terminal_penalty("tree-insertion-terminal-penalty")
    {
      features[feat_penalty] = -1.0;
      features_terminal[feat_terminal_penalty] = -1.0;
      features_terminal[feat_penalty] = -1.0;
      
      attributes["tree-fallback"] = attribute_set_type::int_type(1);
      attributes_terminal["tree-fallback"] = attribute_set_type::int_type(1);
      attributes_terminal["insertion"]     = attribute_set_type::int_type(1);

      if (! goal.empty() || ! non_terminal.empty())
	if (goal.empty() || non_terminal.empty())
	  throw std::runtime_error("You should specify both of goal and non-terminal or none");
    }
    
    TreeGrammarFallback()
      : goal(), non_terminal(),
	feat_penalty("tree-insertion-penalty"),
	feat_terminal_penalty("tree-insertion-terminal-penalty")
    {
      features[feat_penalty] = -1.0;
      features_terminal[feat_terminal_penalty] = -1.0;
      features_terminal[feat_penalty] = -1.0;
      
      attributes["tree-fallback"] = attribute_set_type::int_type(1);
      attributes_terminal["tree-fallback"] = attribute_set_type::int_type(1);
      attributes_terminal["insertion"]     = attribute_set_type::int_type(1);
    }
    
    transducer_ptr_type clone() const { return transducer_ptr_type(new TreeGrammarFallback(*this)); }

    void assign(const hypergraph_type& graph)
    {
      if (non_terminal.empty())
	__assign(graph);
      else
	__assign(graph, goal, non_terminal);
    }
    
  private:
    void __assign(const hypergraph_type& graph, const symbol_type& goal, const symbol_type& non_terminal)
    {
      rules.clear();
      clear();
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = *eiter;
	
	if (! rules.insert(edge.rule).second) continue;
	
	difference_type num_terminal = 0;
	non_terminals.clear();
	symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	for (symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) {
	  non_terminals.push_back(riter->is_non_terminal() ? non_terminal.non_terminal(riter->non_terminal_index()) : *riter);
	  num_terminal += (! riter->is_non_terminal());
	}

	const symbol_type lhs(edge.head == graph.goal ? goal : non_terminal);
	
	rule_ptr_type rule_source(rule_type::create(rule_type(edge.rule->lhs, edge.rule->rhs.begin(), edge.rule->rhs.end())));
	rule_ptr_type rule_target(rule_type::create(rule_type(lhs, non_terminals.begin(), non_terminals.end())));
	
	if (num_terminal) {
	  features_terminal[feat_terminal_penalty] = - num_terminal;
	  
	  insert(rule_pair_type(rule_source, rule_target, features_terminal, attributes_terminal));
	} else
	  insert(rule_pair_type(rule_source, rule_target, features, attributes));
      }
    }
    
    // simply copy and preserve the same non-terminal category...
    void __assign(const hypergraph_type& graph)
    {
      rules.clear();
      clear();
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = *eiter;
	
	if (! rules.insert(edge.rule).second) continue;

	difference_type num_terminal = 0;
	symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	for (symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
	  num_terminal += (! riter->is_non_terminal());
	
	rule_ptr_type rule(rule_type::create(rule_type(edge.rule->lhs, edge.rule->rhs.begin(), edge.rule->rhs.end())));
	
	if (num_terminal) {
	  features_terminal[feat_terminal_penalty] = - num_terminal;
	  
	  insert(rule_pair_type(rule, rule, features_terminal, attributes_terminal));
	} else
	  insert(rule_pair_type(rule, rule, features, attributes));
      }
    }
    
  private:
    symbol_type goal;
    symbol_type non_terminal;
    
    feature_type feat_penalty;
    feature_type feat_terminal_penalty;
    
    graph_rule_ptr_set_type rules;
    feature_set_type        features;
    feature_set_type        features_terminal;
    attribute_set_type      attributes;
    attribute_set_type      attributes_terminal;
    non_terminal_set_type   non_terminals;
  };
  
};

#endif
