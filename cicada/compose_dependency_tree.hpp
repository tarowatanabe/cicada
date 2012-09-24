// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// tranform lattice (actually, sentence) and dependency, into tree structure
//

#ifndef __CICADA__COMPOSE_DEPENDENCY_TREE__HPP__
#define __CICADA__COMPOSE_DEPENDENCY_TREE__HPP__ 1

#include <vector>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/dependency.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/sort_topologically.hpp>

namespace cicada
{
  struct ComposeDependencyTree
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;

    typedef Lattice    lattice_type;
    typedef Dependency dependency_type;
    typedef HyperGraph hypergraph_type;

    typedef std::vector<size_type, std::allocator<size_type> > index_set_type;
    typedef std::vector<index_set_type, std::allocator<index_set_type> > dependency_map_type;
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_map_type;

    ComposeDependencyTree(const symbol_type& __goal,
			  const symbol_type& __non_terminal)
      : goal(__goal), non_terminal(__non_terminal) {}

    void operator()(const lattice_type& lattice,
		    const dependency_type& dependency,
		    hypergraph_type& graph)
    {
      
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tail_set_type;
      typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
      
      if (lattice.size() != dependency.size())
	throw std::runtime_error("invalid lattice and dependency");

      graph.clear();
      
      if (dependency.empty()) return;

      dependency_map.clear();
      dependency_map.resize(dependency.size() + 1);
      
      node_map.clear();
      node_map.resize(dependency.size() + 1, hypergraph_type::invalid);

      tail_set_type tails;
      symbol_set_type rhs;
      
      for (size_type pos = 0; pos != lattice.size(); ++ pos) {
	if (lattice[pos].size() != 1)
	  throw std::runtime_error("this is not a sentential lattice!");

	if (lattice[pos].front().distance != 1)
	  throw std::runtime_error("this is not a sentential lattice!");
	
	dependency_map[dependency[pos]].push_back(pos + 1);
      }
      
      if (dependency_map.front().empty()) return;
      
      tails.clear();
      rhs.clear();

      index_set_type::const_iterator iiter_end = dependency_map.front().end();
      for (index_set_type::const_iterator iiter = dependency_map.front().begin(); iiter != iiter_end; ++ iiter) {
	const size_type antecedent = *iiter;
	
	if (node_map[antecedent] == hypergraph_type::invalid)
	  node_map[antecedent] = graph.add_node().id;
	
	tails.push_back(node_map[antecedent]);
	rhs.push_back(non_terminal);
      }
      
      if (node_map[0] == hypergraph_type::invalid)
	node_map[0] = graph.add_node().id;
      
      hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
      edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(goal, rhs.begin(), rhs.end()));
      
      graph.connect_edge(edge.id, node_map[0]);
      graph.goal = node_map[0];
      
      for (size_type id = 1; id != dependency_map.size(); ++ id) {
	tails.clear();
	rhs.clear();
	
	index_set_type::const_iterator iiter_begin = dependency_map[id].begin();
	index_set_type::const_iterator iiter_end   = dependency_map[id].end();
	index_set_type::const_iterator iiter_lex   = std::lower_bound(iiter_begin, iiter_end, id);
	
	// left...
	for (index_set_type::const_iterator iiter = iiter_begin; iiter != iiter_lex; ++ iiter) {
	  const size_type antecedent = *iiter;
	  
	  if (node_map[antecedent] == hypergraph_type::invalid)
	    node_map[antecedent] = graph.add_node().id;
	  
	  tails.push_back(node_map[antecedent]);
	  rhs.push_back(non_terminal);
	}
	
	// head
	rhs.push_back(lattice[id - 1].front().label);
	
	// right...
	for (index_set_type::const_iterator iiter = iiter_lex; iiter != iiter_end; ++ iiter) {
	  const size_type antecedent = *iiter;
	  
	  if (node_map[antecedent] == hypergraph_type::invalid)
	    node_map[antecedent] = graph.add_node().id;
	  
	  tails.push_back(node_map[antecedent]);
	  rhs.push_back(non_terminal);
	}
	
	if (node_map[id] == hypergraph_type::invalid)
	  node_map[id] = graph.add_node().id;
	
	hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
	edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(non_terminal, rhs.begin(), rhs.end()));
	
	graph.connect_edge(edge.id, node_map[id]);
      }
      
      if (! graph.nodes.empty() && graph.is_valid())
	graph.topologically_sort();
    }
    
    dependency_map_type dependency_map;
    node_map_type node_map;
    
    const symbol_type goal;
    const symbol_type non_terminal;
  };
  
  inline
  void compose_dependency_tree(const Lattice& lattice, const Dependency& dependency, HyperGraph& graph, const Symbol goal="[s]", const Symbol non_terminal="[x]")
  {
    ComposeDependencyTree composer(goal, non_terminal);
    composer(lattice, dependency, graph);
  }

};

#endif
