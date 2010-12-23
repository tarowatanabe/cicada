// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_ALIGNMENT__HPP__
#define __CICADA__COMPOSE_ALIGNMENT__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/span_edge.hpp>

#include <utils/sgi_hash_set.hpp>
#include <utils/sgi_hash_map.hpp>

namespace cicada
{
  struct ComposeAlignment
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;
    
    typedef Lattice    lattice_type;
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;
    
    typedef attribute_set_type::attribute_type attribute_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;
    
    ComposeAlignment(const symbol_type& non_terminal,
		     const grammar_type& __grammar)
      : grammar(__grammar),
	attr_target_word("target-word"),
	attr_source_position("source-position"),
	attr_target_position("target-position")
    {
      // initializer...
      rule_goal = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, non_terminal.non_terminal(1))));
      
      std::vector<symbol_type, std::allocator<symbol_type> > sequence(2);
      sequence.front() = non_terminal.non_terminal(1);
      sequence.back()  = non_terminal.non_terminal(2);
      
      rule_x1_x2 = rule_type::create(rule_type(non_terminal.non_terminal(), sequence.begin(), sequence.end()));
    }

    void operator()(const lattice_type& source, const lattice_type& target, hypergraph_type& graph)
    {
      // here, we assume linear-lattice, mearning sentence!
      graph.clear();
      
      if (source.empty() || target.empty()) return;
      
      hypergraph_type::id_type node_prev = hypergraph_type::invalid;
      
      for (size_t src = 0; src != source.size(); ++ src) {
	hypergraph_type::node_type& node = graph.add_node();
	
	for (int trg = -1; trg < static_cast<int>(target.size()); ++ trg) {
	  const symbol_type& target_symbol = (trg < 0 ? vocab_type::EPSILON : target[trg].front().label);
	  
	  for (size_t id = 0; id != grammar.size(); ++ id) {
	    const transducer_type& transducer = grammar[id];
	    const transducer_type::id_type transducer_node = transducer.next(transducer.root(), target_symbol);
	    
	    if (transducer_node == transducer.root()) continue;
	    
	    const transducer_type::rule_pair_set_type& rules = transducer.rules(transducer_node);
	    
	    if (rules.empty()) continue;
	    
	    transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	    for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	      hypergraph_type::edge_type& edge = graph.add_edge();
	      edge.rule       = riter->target;   // actually, this will be a source word!
	      edge.features   = riter->features;
	      edge.attributes = riter->attributes;
	      
	      edge.attributes[attr_source_position] = attribute_set_type::int_type(src);
	      edge.attributes[attr_target_position] = attribute_set_type::int_type(trg);
	      edge.attributes[attr_target_word] = static_cast<const std::string&>(target_symbol);
	      
	      graph.connect_edge(edge.id, node.id);
	    }
	  }
	}
	
	if (node_prev == hypergraph_type::invalid)
	  node_prev = node.id;
	else {
	  hypergraph_type::id_type tails[2] = {node_prev, node.id};
	  hypergraph_type::edge_type& edge = graph.add_edge(tails, tails + 2);
	  edge.rule = rule_x1_x2;
	  
	  node_prev = graph.add_node().id;
	  graph.connect_edge(edge.id, node_prev);
	}
      }

      hypergraph_type::edge_type& edge = graph.add_edge(&node_prev, (&node_prev) + 1);
      edge.rule = rule_goal;
      
      hypergraph_type::node_type& node = graph.add_node();
      
      graph.connect_edge(edge.id, node.id);
      
      graph.goal = node.id;
      // we do not have to topologically sort this!
    }
    
    void operator()(const hypergraph_type& source, const lattice_type& target, hypergraph_type& graph)
    {
      typedef std::pair<int, int> span_type;
      typedef std::vector<span_type, std::allocator<span_type> > span_set_type;

      graph.clear();

      if (! source.is_valid() || target.empty()) return;
      
      // we will simply construct graph sharing the same as source,
      // but differ in that we have extra edges...
      // we assume penn-treebank style tree...

      span_set_type spans(source.edges.size());
      cicada::span_edge(source, spans);

      graph.goal = source.goal;
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	hypergraph_type::node_type& node = graph.add_node();

	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = niter->edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = niter->edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  
	  if (edge_source.rule->rhs.size() != 1 || edge_source.rule->rhs.front().is_non_terminal()) {
	    // simply copy...
	    hypergraph_type::edge_type& edge = graph.add_edge(edge_source.tails.begin(), edge_source.tails.end());
	    edge.rule       = edge_source.rule;
	    edge.features   = edge_source.features;
	    edge.attributes = edge_source.attributes;
	    
	    graph.connect_edge(edge.id, node.id);
	  } else {
	    // we are terminal... derive source-positon by span
	    const int src = spans[edge_source.id].first;
	    
	    for (int trg = -1; trg < static_cast<int>(target.size()); ++ trg) {
	      const symbol_type& target_symbol = (trg < 0 ? vocab_type::EPSILON : target[trg].front().label);
	      
	      for (size_t id = 0; id != grammar.size(); ++ id) {
		const transducer_type& transducer = grammar[id];
		const transducer_type::id_type transducer_node = transducer.next(transducer.root(), target_symbol);
		
		if (transducer_node == transducer.root()) continue;
		
		const transducer_type::rule_pair_set_type& rules = transducer.rules(transducer_node);
		
		if (rules.empty()) continue;
		
		transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
		for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
		  hypergraph_type::edge_type& edge = graph.add_edge();
		  edge.rule       = rule_type::create(rule_type(edge_source.rule->lhs, riter->target->rhs));
		  edge.features   = riter->features + edge_source.features;
		  edge.attributes = riter->attributes + edge_source.attributes;
		  
		  // we need to compute source-pos!
		  edge.attributes[attr_source_position] = attribute_set_type::int_type(src);
		  edge.attributes[attr_target_position] = attribute_set_type::int_type(trg);
		  edge.attributes[attr_target_word] = static_cast<const std::string&>(target_symbol);
		  
		  graph.connect_edge(edge.id, node.id);
		}
	      }
	    }
	  }
	}
      }
    }
    
    const grammar_type& grammar;
    
    const attribute_type attr_target_word;
    const attribute_type attr_source_position;
    const attribute_type attr_target_position;
    
    rule_ptr_type rule_goal;
    rule_ptr_type rule_x1_x2;
  };
  
};

#endif
