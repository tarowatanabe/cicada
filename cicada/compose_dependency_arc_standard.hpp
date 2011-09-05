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
#include <cicada/hypergraph.hpp>
#include <cicada/sort_topologically.hpp>
#include <cicada/remove_epsilon.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/bithack.hpp>
#include <utils/chart.hpp>

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
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef attribute_set_type::attribute_type attribute_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    ComposeDependencyArcStandard()
      : attr_dependency_pos("dependency-pos"),
	attr_dependency_head("dependency-head"),
	attr_dependency_dependent("dependency-dependent")
    {
      node_map.set_empty_key(id_type(-1));

      rule_goal   = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, vocab_type::X)));
      rule_reduce = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(2, vocab_type::X)));
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
    
    typedef google::dense_hash_map<id_type, hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<id_type> > node_map_type;
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      actives.clear();
      actives.resize(lattice.size() + 1);
      
      // initialize actives by axioms... (terminals)
      id_type id = 1;
      for (size_t pos = 0; pos != lattice.size(); ++ pos) {
	lattice_type::arc_set_type::const_iterator aiter_end = lattice[pos].end();
	for (lattice_type::arc_set_type::const_iterator aiter = lattice[pos].begin(); aiter != aiter_end; ++ aiter, ++ id) {
	  hypergraph_type::edge_type& edge = graph.add_edge();
	  edge.rule = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(1, aiter->label)));
	  edge.features = aiter->features;
	  edge.attributes[attr_dependency_pos] = attribute_set_type::int_type(id);
	  
	  const hypergraph_type::id_type node_id = graph.add_node().id;
	  
	  graph.connect_edge(edge.id, node_id);
	  
	  actives(pos, pos + aiter->distance).push_back(item_type(id, node_id));
	}
      }
      
      hypergraph_type::edge_type::node_set_type tails(2);
      
      for (size_t length = 2; length <= lattice.size(); ++ length)
	for (size_t first = 0; first + length <= lattice.size(); ++ first) {
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
		
		{
		  // left attachment
		  hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		  edge.rule = rule_reduce;
		  edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(riter->id);
		  edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(liter->id);
		  
		  std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(riter->id, 0));
		  if (result.second) {
		    result.first->second = graph.add_node().id;
		    cell.push_back(item_type(riter->id, result.first->second));
		  }
		  
		  graph.connect_edge(edge.id, result.first->second);
		}
		
		{
		  // right attachment
		  hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		  edge.rule = rule_reduce;
		  edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(liter->id);
		  edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(riter->id);
		  
		  std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(liter->id, 0));
		  if (result.second) {
		    result.first->second = graph.add_node().id;
		    cell.push_back(item_type(liter->id, result.first->second));
		  }
		  
		  graph.connect_edge(edge.id, result.first->second);
		}
	      }
	  }
	}
      
      // goal
      
      const hypergraph_type::id_type goal_id = graph.add_node().id;
      graph.goal = goal_id;
      
      const item_set_type& goals = actives(0, lattice.size());
      
      item_set_type::const_iterator giter_end = goals.end();
      for (item_set_type::const_iterator giter = goals.begin(); giter != giter_end; ++ giter) {
	hypergraph_type::edge_type& edge = graph.add_edge(&(giter->node), &(giter->node) + 1);
	edge.rule = rule_goal;
	edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(0);
	edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(giter->id);
	
	graph.connect_edge(edge.id, goal_id);
      }
      
      graph.topologically_sort();
    }
    
  private:
    const attribute_type attr_dependency_pos;
    const attribute_type attr_dependency_head;
    const attribute_type attr_dependency_dependent;

    active_chart_type     actives;
    node_map_type         node_map;
    
    rule_ptr_type rule_goal;
    rule_ptr_type rule_reduce;
  };

  inline
  void compose_dependency_arc_standard(const Lattice& lattice, HyperGraph& graph)
  {
    ComposeDependencyArcStandard composer;
    composer(lattice, graph);
  }
};

#endif
