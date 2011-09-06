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
	  
	  // is this correct...?
	  // we need to shift + 1 for correct indexing...
	  // [h3, j, h3 j, j + 1] where j == pos + 1 and h3 should starts from -1
	  
	  const int j = pos + 1;
	  for (int h3 = -1; h3 != j; ++ h3) {
	    //std::cerr << "axiom [" << h3 << ", " << j << ", " << h3 << " " << j << ", " << j + aiter->distance << ']' << std::endl;
	    actives(j, j + aiter->distance)(h3 + 1, h3 + 1, j + 1) = node_id;
	  }
	}
      }
      
      hypergraph_type::edge_type::node_set_type tails(2);
      
      const int last_max = lattice.size() + 1;
      for (int last = 2; last <= last_max; ++ last) 
	for (int length = 2; last - length >= 0; ++ length) {
	  const int first = last - length;

	  item_set_type& items = actives(first, last);
	  
	  //
	  // is it correct???
	  // [h1, first, h2 h3, middle] [h3, middle, h4 h5, last]
	  //
	  // first  <= h3 < middle
	  // middle <= h5 < last
	  // h3 <= h4 < h5
	  // -1 <= h1 < h3 (or first??? given h3 < middle...?)
	  // h1 <= h2 < h3
	  
	  for (int middle = first + 1; middle < last; ++ middle) {
	    const item_set_type& items_left  = actives(first, middle);
	    const item_set_type& items_right = actives(middle, last);

	    //std::cerr << "span: " << first << ".." << middle << " " << middle << ".." << last << std::endl;
	    
	    for (int h3 = first; h3 < middle; ++ h3)
	      for (int h5 = middle; h5 < last; ++ h5)
		for (int h4 = h3; h4 < h5; ++ h4)
		  for (int h1 = -1; h1 < first; ++ h1)
		    for (int h2 = h1; h2 < h3; ++ h2) {
		      const hypergraph_type::id_type item_left  = items_left(h1 + 1, h2 + 1, h3 + 1);
		      const hypergraph_type::id_type item_right = items_right(h3 + 1, h4 + 1, h5 + 1);
		      
		      const bool item_left_epsilon = (first == 0 && middle == 1 && h1 == -1 && h2 == -1 && h3 == 0);
		      const bool item_left_valid  = item_left_epsilon || item_left != hypergraph_type::invalid;
		      const bool item_right_valid = item_right != hypergraph_type::invalid;
		      
		      if (! item_left_valid || ! item_right_valid) continue;
		      
		      //std::cerr << "item1 [" << h1 << ", " << first << ", " << h2 << " " << h3 << ", " << middle << ']' << std::endl;
		      //std::cerr << "item2 [" << h3 << ", " << middle << ", " << h4 << " " << h5 << ", " << last << ']' << std::endl;
		      
		      tails.front() = item_left;
		      tails.back()  = item_right;
		      
		      // [h1, i, h2 h5, j] (la1; h5 -> h4)
		      if (h4 > 0) {
			// left attachment
			hypergraph_type::edge_type& edge = graph.add_edge(tails.begin() + item_left_epsilon, tails.end());
			edge.rule = (item_left_epsilon ? rule_reduce1 : rule_reduce2);
			edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(h5);
			edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(h4);
			
			hypergraph_type::id_type& cell =  items(h1 + 1, h2 + 1, h5 + 1);
			if (cell == hypergraph_type::invalid)
			  cell = graph.add_node().id;
			
			graph.connect_edge(edge.id, cell);
		      }
		      
		      // [h1, i, h2 h4, j] (ra1; h4 -> h5)
		      if (h4 >= 0) {
			// right attachment
			hypergraph_type::edge_type& edge = graph.add_edge(tails.begin() + item_left_epsilon, tails.end());
			edge.rule = (item_left_epsilon ? rule_reduce1 : rule_reduce2);
			edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(h4);
			edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(h5);
			
			hypergraph_type::id_type& cell =  items(h1 + 1, h2 + 1, h4 + 1);
			if (cell == hypergraph_type::invalid)
			  cell = graph.add_node().id;
			
			graph.connect_edge(edge.id, cell);
		      } 
		      
		      // [h1, i, h4 h5, j] (la2; h5 -> h2)
		      if (h2 > 0) {
			// left attachment
			hypergraph_type::edge_type& edge = graph.add_edge(tails.begin() + item_left_epsilon, tails.end());
			edge.rule = (item_left_epsilon ? rule_reduce1 : rule_reduce2);
			edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(h5);
			edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(h2);
			
			hypergraph_type::id_type& cell =  items(h1 + 1, h4 + 1, h5 + 1);
			if (cell == hypergraph_type::invalid)
			  cell = graph.add_node().id;
			
			graph.connect_edge(edge.id, cell);
		      }
		      
		      // [h1, i, h2 h4, j] (ra2; h2 -> h5)
		      if (h2 >= 0) {
			// right attachment
			hypergraph_type::edge_type& edge = graph.add_edge(tails.begin() + item_left_epsilon, tails.end());
			edge.rule = (item_left_epsilon ? rule_reduce1 : rule_reduce2);
			edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(h2);
			edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(h5);
			
			hypergraph_type::id_type& cell =  items(h1 + 1, h2 + 1, h4 + 1);
			if (cell == hypergraph_type::invalid)
			  cell = graph.add_node().id;
			
			graph.connect_edge(edge.id, cell);
		      }
		    }
	  }
	}
      
      // final...
      graph.goal = actives(0, last_max)(-1 + 1, -1 + 1, 0 + 1);
      if (graph.is_valid())
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
