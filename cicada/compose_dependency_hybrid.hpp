// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_DEPENDENCY_HYBRID__HPP__
#define __CICADA__COMPOSE_DEPENDENCY_HYBRID__HPP__ 1

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

namespace cicada
{
  // hybrid parser based on the deduction system presented in
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
  // }
  //
  // which is originally presented in
  //
  // @INPROCEEDINGS{Yamada03statisticaldependency,
  //   author = {Hiroyasu Yamada and Yuji Matsumoto},
  //   title = {Statistical Dependency Analysis with Support Vector Machines},
  //   booktitle = {In Proceedings of IWPT},
  //   year = {2003},
  //   pages = {195--206}
  // }

  struct ComposeDependencyHybrid
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

    ComposeDependencyHybrid()
      : attr_dependency_pos("dependency-pos"),
	attr_dependency_head("dependency-head"),
	attr_dependency_dependent("dependency-dependent")
    {
      rule_epsilon = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(1, vocab_type::EPSILON)));
      rule_x1_x2   = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(2, vocab_type::X)));
    }
    
    typedef utils::chart<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> >  active_chart_type;
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      actives.clear();
      actives.resize(lattice.size() + 2, hypergraph_type::invalid);
      
      // initialize actives by axioms... (terminals)
      
      // root...
      // we will insert pseudo edge, but this will be "removed"
      
      hypergraph_type::edge_type& edge = graph.add_edge();
      edge.rule = rule_epsilon;
      edge.attributes[attr_dependency_pos] = attribute_set_type::int_type(0);
      
      const hypergraph_type::id_type node_id = graph.add_node().id;
      
      graph.connect_edge(edge.id, node_id);
      
      actives(0, 1) = node_id;

      if (edge.id != 0)
	throw std::runtime_error("invalid edge id?");
      if (node_id != 0)
	throw std::runtime_error("invalid node id?");
      
      // we assume sentence-like data, not general lattice!
      
      for (size_t pos = 0; pos != lattice.size(); ++ pos) {
	
	if (lattice[pos].size() != 1)
	  throw std::runtime_error("this is not a sentential lattice!");
	
	// here, we will construct a partial hypergraph...
	lattice_type::arc_set_type::const_iterator aiter_end  = lattice[pos].end();
	for (lattice_type::arc_set_type::const_iterator aiter = lattice[pos].begin(); aiter != aiter_end; ++ aiter) {
	  hypergraph_type::edge_type& edge = graph.add_edge();
	  edge.rule = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(1, aiter->label)));
	  
	  edge.features = aiter->features;
	  edge.attributes[attr_dependency_pos] = attribute_set_type::int_type(pos + 1);
	  
	  const hypergraph_type::id_type node_id = graph.add_node().id;
	  
	  graph.connect_edge(edge.id, node_id);
	  
	  actives(pos + 1, pos + aiter->distance + 1) = node_id;
	}
      }
      
      hypergraph_type::edge_type::node_set_type tails(2);
      
      for (int last = 2; last <= static_cast<int>(lattice.size() + 1); ++ last) 
	for (int length = 2; last - length >= 0; ++ length) {
	  const int first = last - length;
	  
	  hypergraph_type::id_type& cell = actives(first, last);
	  
	  for (int middle = first + 1; middle < last; ++ middle) 
	    if (actives(first, middle) != hypergraph_type::invalid && actives(middle, last) != hypergraph_type::invalid) {
	      
	      if (cell == hypergraph_type::invalid)
		cell = graph.add_node().id;
	      
	      tails.front() = actives(first, middle);
	      tails.back()  = actives(middle, last);
	      
	      if (last < static_cast<int>(lattice.size() + 1)) {
		// left attachment
		hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		edge.rule = rule_x1_x2;
		edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(last);
		edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(middle);
		
		graph.connect_edge(edge.id, cell);
	      }

	      
	      {
		// right attachment
		hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		edge.rule = rule_x1_x2;
		edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(first);
		edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(middle);
		
		graph.connect_edge(edge.id, cell);
	      }
	    }
	}
      
      hypergraph_type::id_type goal_id = actives(0, lattice.size() + 1);
      
      if (goal_id == hypergraph_type::invalid) return;
      
      graph.goal = goal_id;
      
      cicada::remove_epsilon(graph);
    }
    
  private:
    const attribute_type attr_dependency_pos;
    const attribute_type attr_dependency_head;
    const attribute_type attr_dependency_dependent;
    
    // we need to keep track of two actives in the first and the second
    active_chart_type     actives;
    
    rule_ptr_type rule_epsilon;
    rule_ptr_type rule_x1_x2;
  };

  inline
  void compose_dependency_hybrid(const Lattice& lattice, HyperGraph& graph)
  {
    ComposeDependencyHybrid composer;
    composer(lattice, graph);
  }
};

#endif
