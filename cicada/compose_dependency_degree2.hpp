// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_DEPENDENCY_DEGREE2__HPP__
#define __CICADA__COMPOSE_DEPENDENCY_DEGREE2__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/sort_topologically.hpp>
#include <cicada/remove_epsilon.hpp>

#include <utils/bithack.hpp>
#include <utils/chart.hpp>
#include <utils/vector3.hpp>

namespace cicada
{
  // degree-2 non-projective parser based on the deduction system presented in
  //
  // @InProceedings{cohen-gomezrodriguez-satta:2011:EMNLP,
  //   author    = {Cohen, Shay B.  and  G\'{o}mez-Rodr\'{i}guez, Carlos  and  Satta, Giorgio},
  //   title     = {Exact Inference for Generative Probabilistic Non-Projective Dependency Parsing},
  //   booktitle = {Proceedings of the 2011 Conference on Empirical Methods in Natural Language Processing},
  //   month     = {July},
  //   year      = {2011},
  //   address   = {Edinburgh, Scotland, UK.},
  //   publisher = {Association for Computational Linguistics},
  //   pages     = {1234--1245},
  //   url       = {http://www.aclweb.org/anthology/D11-1114}
  // }
  //

  struct ComposeDependencyDegree2
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

    ComposeDependencyDegree2()
      : attr_dependency_pos("dependency-pos"),
	attr_dependency_head("dependency-head"),
	attr_dependency_dependent("dependency-dependent")
    {
      rule_reduce1 = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(1, vocab_type::X)));
      rule_reduce2 = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(2, vocab_type::X)));
    }
    
    // for now, we use chart structure + array3, but I think we can easily optimize it away given
    // the diagonal structure employed in the tabulation
    //
    // [h_1, i, h_2 h_3, j]
    // where i <= h_2 < h_3 < j
    // h_1 < i (??) (h_1 can be -1..????)
    //
    // since h_1 can be -1, we need to shift its index... (?)
    
    typedef utils::vector3<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > item_set_type;
    typedef utils::chart<item_set_type, std::allocator<item_set_type> >  active_chart_type;
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      actives.clear();
      actives.resize(lattice.size() + 2, item_set_type(lattice.size() + 2, lattice.size() + 2, lattice.size() + 2, hypergraph_type::invalid));
      
      // initialize actives by axioms... (terminals)
      
      for (size_t pos = 0; pos != lattice.size(); ++ pos) {
	
	if (lattice[pos].size() != 1)
	  throw std::runtime_error("this is not a sentential lattice!");
	
	// here, we will construct a partial hypergraph...
	lattice_type::arc_set_type::const_iterator aiter_end  = lattice[pos].end();
	for (lattice_type::arc_set_type::const_iterator aiter = lattice[pos].begin(); aiter != aiter_end; ++ aiter) {
	  
	  if (aiter->distance != 1)
	    throw std::runtime_error("this is not a sentential lattice");
	  
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
      
      const int last_max = lattice.size() + 1;
      for (int last = 2; last <= last_max; ++ last) 
	for (int length = 2; last - length >= 0; ++ length) {
	  const int first = last - length;
	  
	  hypergraph_type::id_type& cell = actives(first, last);
	  
	  cell = graph.add_node().id;
	  
	  for (int middle = first + 1; middle < last; ++ middle) {
	    tails.front() = actives(first, middle);
	    tails.back()  = actives(middle, last);
	    
	    if (first == 0 && middle == 1) {
	      if (last < last_max) {
		hypergraph_type::edge_type& edge = graph.add_edge(tails.begin() + 1, tails.end());
		edge.rule = rule_reduce1;
		edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(last);
		edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(middle);
		
		graph.connect_edge(edge.id, cell);
	      }
	      
	      {
		hypergraph_type::edge_type& edge = graph.add_edge(tails.begin() + 1, tails.end());
		edge.rule = rule_reduce1;
		edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(first);
		edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(middle);
		
		graph.connect_edge(edge.id, cell);
	      }
	    } else {
	      if (last < last_max) {
		// left attachment
		hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		edge.rule = rule_reduce2;
		edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(last);
		edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(middle);
		
		graph.connect_edge(edge.id, cell);
	      }
	      
	      {
		// right attachment
		hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		edge.rule = rule_reduce2;
		edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(first);
		edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(middle);
		
		graph.connect_edge(edge.id, cell);
	      }
	    }
	  }
	}
      
      // final...
      graph.goal = actives(0, last_max);
      graph.topologically_sort();
    }
    
  private:
    const attribute_type attr_dependency_pos;
    const attribute_type attr_dependency_head;
    const attribute_type attr_dependency_dependent;
    
    // we need to keep track of two actives in the first and the second
    active_chart_type     actives;
    
    rule_ptr_type rule_reduce1;
    rule_ptr_type rule_reduce2;
  };

  inline
  void compose_dependency_degree2(const Lattice& lattice, HyperGraph& graph)
  {
    ComposeDependencyDegree2 composer;
    composer(lattice, graph);
  }
};

#endif
