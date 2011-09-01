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
#include <utils/sgi_hash_set.hpp>

#include <google/dense_hash_map>

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
      node_map.set_empty_key(id_type(-1));
      
      rule_epsilon = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(1, vocab_type::EPSILON)));
      rule_x1_x2   = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(2, vocab_type::X)));
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

#ifdef HAVE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_set<size_type, utils::hashmurmur<size_t>, std::equal_to<size_type>, std::allocator<size_type> > pos_set_type;
#else
    typedef sgi::hash_set<size_type, utils::hashmurmur<size_t>, std::equal_to<size_type>, std::allocator<size_type> > pos_set_type;
#endif
    typedef std::vector<pos_set_type, std::allocator<pos_set_type> > pos_map_type;
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      actives.clear();
      actives.resize(lattice.size() + 2);
      
      // initialize actives by axioms... (terminals)
      
      // root...
      // we will insert pseudo edge, but this will be "removed"
      
      hypergraph_type::edge_type& edge = graph.add_edge();
      edge.rule = rule_epsilon;
      edge.attributes[attr_dependency_pos] = attribute_set_type::int_type(0);
      
      const hypergraph_type::id_type node_id = graph.add_node().id;
      
      graph.connect_edge(edge.id, node_id);
      
      actives_first(0, 1).push_back(item_type(0, node_id));

      if (edge.id != 0)
	throw std::runtime_error("invalid edge id?");
      if (node_id != 0)
	throw std::runtime_error("invalid node id?");
      
      // first, compute inverse index...
      
      pos_map_type pos_map(lattice.size() + 1);
      for (size_t pos = 0; pos != lattice.size(); ++ pos) {
	lattice_type::arc_set_type::const_iterator aiter_end  = lattice[pos].end();
	for (lattice_type::arc_set_type::const_iterator aiter = lattice[pos].begin(); aiter != aiter_end; ++ aiter)
	  pos_map[pos + aiter->distance].insert(pos);
      }
      
      id_type id = 1;
      for (size_t pos = 0; pos != lattice.size(); ++ pos) {
	// here, we will construct a partial hypergraph...
	lattice_type::arc_set_type::const_iterator aiter_end  = lattice[pos].end();
	for (lattice_type::arc_set_type::const_iterator aiter = lattice[pos].begin(); aiter != aiter_end; ++ aiter, ++ id) {
	  
	  hypergraph_type::edge_type& edge = graph.add_edge();
	  edge.rule = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(1, aiter->label)));
	  
	  edge.features = aiter->features;
	  edge.attributes[attr_dependency_pos] = attribute_set_type::int_type(id);
	  
	  const hypergraph_type::id_type node_id = graph.add_node().id;
	  
	  graph.connect_edge(edge.id, node_id);
	  
	  actives_first(pos + 1, pos + aiter->distance + 1).push_back(item_type(id, node_id));
	  
	  pos_set_type::const_iterator piter_end = pos_map[pos].end();
	  for (pos_set_type::const_iterator piter = pos_map[pos].end(); piter != piter_end; ++ piter)
	    actives_second(*piter + 1, pos + 1).push_back(item_type(id, node_id));
	}
      }

      hypergraph_type::edge_type::node_set_type tails(2);
      rule_type::symbol_set_type                rhs(2);

      for (int last = 2; last <= static_cast<int>(lattice.size() + 1); ++ last) 
	for (int length = 2; last - length >= 0; ++ length) {
	  const int first = last - length;
	  
	  node_map.clear();
	  
	  item_set_type& cell_first  = actives_first(first, last);
	  item_set_type& cell_second = actives_second(first, last);
	  
	  for (int middle = first + 1; middle < last; ++ middle) {
	    	    
	    {
	      // second: left-attachment
	      const item_set_type& items_left  = actives_second(first, middle);
	      const item_set_type& items_right = actives_second(middle, last);
	      
	      if (! items_left.empty() && ! items_right.empty()) {
		item_set_type::const_iterator liter_begin = items_left.begin();
		item_set_type::const_iterator liter_end   = items_left.end();
		item_set_type::const_iterator riter_begin = items_right.begin();
		item_set_type::const_iterator riter_end   = items_right.end();
		
		for (item_set_type::const_iterator liter = liter_begin; liter != liter_end; ++ liter)
		  for (item_set_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter) {
		    tails.front() = liter->node;
		    tails.back()  = riter->node;
		    
		    const symbol_type& lhs = rhs.back();
		    
		    hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		    edge.rule = rule_type::create(rule_type(lhs, rhs));
		    edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(riter->id);
		    edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(liter->id);
		    
		    std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(lhs, riter->id), 0));
		    if (result.second) {
		      result.first->second = graph.add_node().id;
		      
		      cell_second.push_back(item_type(riter->id, result.first->second));
		    }
		    
		    graph.connect_edge(edge.id, result.first->second);
		  }
	      }
	    }
	    
	    
	    {
	      // first: right-attachment
	      const item_set_type& items_left  = actives_first(first, middle);
	      const item_set_type& items_right = actives_first(middle, last);
	      
	      if (! items_left.empty() && ! items_right.empty()) {
		item_set_type::const_iterator liter_begin = items_left.begin();
		item_set_type::const_iterator liter_end   = items_left.end();
		item_set_type::const_iterator riter_begin = items_right.begin();
		item_set_type::const_iterator riter_end   = items_right.end();
		
		for (item_set_type::const_iterator liter = liter_begin; liter != liter_end; ++ liter)
		  for (item_set_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter) {
		    tails.front() = liter->node;
		    tails.back()  = riter->node;
		    
		    const symbol_type& lhs = rhs.front();
		    
		    hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		    edge.rule = rule_type::create(rule_type(lhs, rhs));
		    edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(liter->id);
		    edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(riter->id);
		    
		    std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(lhs, liter->id), 0));
		    if (result.second) {
		      result.first->second = graph.add_node().id;
		      
		      cell_first.push_back(item_type(liter->id, result.first->second));
		    }
		    
		    graph.connect_edge(edge.id, result.first->second);
		  }
	      }
	    }
	  }
	}
      
      
      // add goals!
      const item_set_type& goals = actives_first(0, lattice.size() + 1);
      
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
    const attribute_type attr_dependency_pos;
    const attribute_type attr_dependency_head;
    const attribute_type attr_dependency_dependent;
    
    // we need to keep track of two actives in the first and the second
    active_chart_type     actives;
    node_map_type         node_map;
    
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
