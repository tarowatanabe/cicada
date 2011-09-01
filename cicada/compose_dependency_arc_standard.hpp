// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_DEPENDENCY_ARC_STANDARD__HPP__
#define __CICADA__COMPOSE_DEPENDENCY_ARC_STANDARD__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/sort_topologically.hpp>
#include <cicada/remove_epsilon.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/bithack.hpp>

#include <google/dense_hash_map>

namespace cicada
{
  // arc standard parser based on the deduction system presented in
  //
  // @InProceedings{kuhlmann-gomezrodriguez-satta:2011:ACL-HLT2011,
  //   author    = {Kuhlmann, Marco  and  G\'{o}mez-Rodr\'{i}guez, Carlos  and  Satta, Giorgio},
  //   title     = {Dynamic Programming Algorithms for Transition-Based Dependency Parsers},
  //   booktitle = {Proceedings of the 49th Annual Meeting of the Association for Computational Linguistics: Human Language Technologies},
  //   month     = {June},
  //   year      = {2011},
  //   address   = {Portland, Oregon, USA},
  //   publisher = {Association for Computational Linguistics},
  //   pages     = {673--682},
  //   url       = {http://www.aclweb.org/anthology/P11-1068}
  //  }
  //

  struct ComposeDependencyArcStandard
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

    ComposeDependencyArcStandard(const grammar_type& __grammar,
				 const bool __pos_mode=false)
      : grammar(__grammar), pos_mode(__pos_mode),
	attr_dependency_pos("dependency-pos"),
	attr_dependency_head("dependency-head"),
	attr_dependency_dependent("dependency-dependent")
    {
      node_map.set_empty_key(symbol_id_type());
    }

    typedef uint32_t id_type;
    
    struct Item
    {
      Item() : id(0), node(hypergraph_type::invalid) {}
      Item(const id_type& __id, const hypergraph_type::id_type& __node)
	: id(__id), node(__node) {}
      
      id_type id;
      hypergraph_type::id_type node;
    };
    
    typedef Item item_type;
    typedef utils::chunk_vector<item_type, 4096 / sizeof(item_type), std::allocator<item_type> > item_set_type;
    typedef utils::chart<item_set_type, std::allocator<item_set_type> >  active_chart_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
    
    typedef std::pair<symbol_type, id_type> symbol_id_type;
    typedef google::dense_hash_map<symbol_id_type, hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<symbol_id_type> > node_map_type;
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      actives.clear();
      non_terminals.clear();
      
      actives.resize(lattice.size() + 2);
      
      // initialize actives by axioms... (terminals)
      
      // root...
      // we will insert pseudo edge, but this will be "removed"
      
      hypergraph_type::edge_type& edge = graph.add_edge();
      edge.rule = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, vocab_type::EPSILON)));
      
      edge.attributes[attr_dependency_pos] = attribute_set_type::int_type(0);
      
      const hypergraph_type::id_type node_id = graph.add_node().id;
      non_terminals.push_back(vocab_type::GOAL);
      
      graph.connect_edge(edge.id, node_id);
      
      actives(0, 1).push_back(item_type(0, node_id));

      if (edge.id != 0)
	throw std::runtime_error("invalid edge id?");
      if (node_id != 0)
	throw std::runtime_error("invalid node id?");

      id_type id = 1;
      for (size_t pos = 0; pos != lattice.size(); ++ pos) {
	// here, we will construct a partial hypergraph...
	
	if (pos_mode) {
	  lattice_type::arc_set_type::const_iterator aiter_end = lattice[pos].end();
	  for (lattice_type::arc_set_type::const_iterator aiter = lattice[pos].begin(); aiter != aiter_end; ++ aiter, ++ id) {
	    const symbol_type terminal = aiter->label.terminal();
	    
	    symbol_type tag = aiter->label.pos();
	    if (tag.empty())
	      tag = vocab_type::X;
	    
	    hypergraph_type::edge_type& edge = graph.add_edge();
	    edge.rule = rule_type::create(rule_type(tag, rule_type::symbol_set_type(1, terminal)));
	    
	    edge.features = aiter->features;
	    edge.attributes[attr_dependency_pos] = attribute_set_type::int_type(id);
	    
	    const hypergraph_type::id_type node_id = graph.add_node().id;
	    non_terminals.push_back(tag);
	    
	    graph.connect_edge(edge.id, node_id);
	    
	    actives(pos + 1, pos + aiter->distance + 1).push_back(item_type(id, node_id));
	  }
	} else {
	  node_map.clear();
	  
	  lattice_type::arc_set_type::const_iterator aiter_end = lattice[pos].end();
	  for (lattice_type::arc_set_type::const_iterator aiter = lattice[pos].begin(); aiter != aiter_end; ++ aiter, ++ id) {
	    // enumerate grammar...
	    
	    for (size_t table = 0; table != grammar.size(); ++ table) {
	      const transducer_type& transducer = grammar[table];
	      
	      const transducer_type::id_type node = transducer.next(transducer.root(), aiter->label);
	      if (node == transducer.root()) continue;
	      
	      const transducer_type::rule_pair_set_type& rules = transducer.rules(node);
	      
	      if (rules.empty()) continue;
	      
	      transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	      for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
		const symbol_type& lhs = riter->source->lhs;
		
		std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(lhs, id), 0));
		if (result.second) {
		  result.first->second = graph.add_node().id;
		  non_terminals.push_back(lhs);
		  
		  actives(pos + 1, pos + aiter->distance + 1).push_back(item_type(id, result.first->second));
		}
		
		hypergraph_type::edge_type& edge = graph.add_edge();
		edge.rule = rule_type::create(rule_type(lhs, rule_type::symbol_set_type(1, aiter->label)));
		
		edge.features = riter->features + aiter->features;
		edge.attributes = riter->attributes;
		
		edge.attributes[attr_dependency_pos] = attribute_set_type::int_type(id);
		
		graph.connect_edge(edge.id, result.first->second);
	      }
	    }
	  }
	}
      }

      hypergraph_type::edge_type::node_set_type tails(2);
      rule_type::symbol_set_type                rhs(2);
      
      for (size_t length = 2; length <= lattice.size() + 1; ++ length)
	for (size_t first = 0; first + length <= lattice.size() + 1; ++ first) {
	  const size_t last = first + length;

	  node_map.clear();
	  item_set_type& cell = actives(first, last);
	  
	  for (size_t middle = first + 1; middle < last; ++ middle) {
	    const item_set_type& items_left  = actives(first, middle);
	    const item_set_type& items_right = actives(middle, last);
	    
	    if (items_left.empty() || items_right.empty()) continue;

	    item_set_type::const_iterator liter_begin = items_left.begin();
	    item_set_type::const_iterator liter_end   = items_left.end();
	    item_set_type::const_iterator riter_begin = items_right.begin();
	    item_set_type::const_iterator riter_end   = items_right.end();
	    
	    for (item_set_type::const_iterator liter = liter_begin; liter != liter_end; ++ liter)
	      for (item_set_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter) {
		tails.front() = liter->node;
		tails.back()  = riter->node;
		
		rhs.front() = non_terminals[liter->node];
		rhs.back()  = non_terminals[riter->node];
		
		if (liter->id) {
		  // left attachment
		  const symbol_type& lhs = rhs.back();

		  hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		  edge.rule = rule_type::create(rule_type(lhs, rhs));
		  edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(riter->id);
		  edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(liter->id);
		  
		  std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(lhs, riter->id), 0));
		  if (result.second) {
		    result.first->second = graph.add_node().id;
		    non_terminals.push_back(lhs);
		    
		    cell.push_back(item_type(riter->id, result.first->second));
		  }

		  graph.connect_edge(edge.id, result.first->second);
		}
		
		if (riter->id) {
		  // right attachment
		  const symbol_type& lhs = rhs.front();

		  hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		  edge.rule = rule_type::create(rule_type(lhs, rhs));
		  edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(liter->id);
		  edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(riter->id);
		  
		  std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(lhs, liter->id), 0));
		  if (result.second) {
		    result.first->second = graph.add_node().id;
		    non_terminals.push_back(lhs);
		    
		    cell.push_back(item_type(liter->id, result.first->second));
		  }
		  
		  graph.connect_edge(edge.id, result.first->second);
		}
	      }
	  }
	}
      
      // add goals!
      const item_set_type& goals = actives(0, lattice.size() + 1);
      
      hypergraph_type::id_type goal_id = hypergraph_type::invalid;
      size_t num_goal = 0;
      item_set_type::const_iterator giter_end = goals.end();
      for (item_set_type::const_iterator giter = goals.begin(); giter != giter_end; ++ giter) {
	num_goal += (giter->id == 0);
	goal_id = utils::bithack::branch(giter->id == 0, giter->node, goal_id);
      }
      
      if (num_goal == 0) return;
      if (num_goal > 1)
	throw std::runtime_error("invalid dependency forest?");
      
      graph.goal = goal_id;
      
      cicada::remove_epsilon(graph);
    }
    
  private:
    const grammar_type& grammar;
    const bool pos_mode;

    const attribute_type attr_dependency_pos;
    const attribute_type attr_dependency_head;
    const attribute_type attr_dependency_dependent;
    
    active_chart_type     actives;
    non_terminal_set_type non_terminals;
    node_map_type         node_map;
  };

  inline
  void compose_dependency_arc_standard(const Grammar& grammar, const Lattice& lattice, HyperGraph& graph, const bool pos_mode=false)
  {
    ComposeDependencyArcStandard(grammar, pos_mode)(lattice, graph);
  }
};

#endif
