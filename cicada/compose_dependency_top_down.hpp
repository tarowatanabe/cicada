// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_DEPENDENCY_TOP_DOWN__HPP__
#define __CICADA__COMPOSE_DEPENDENCY_TOP_DOWN__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/sort_topologically.hpp>
#include <cicada/remove_epsilon.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/bithack.hpp>
#include <utils/bit_vector.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/sgi_hash_map.hpp>

#include <boost/fusion/tuple.hpp>

namespace cicada
{
  // top-down parser 

  struct ComposeDependencyTopDown
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

    ComposeDependencyTopDown()
      : attr_dependency_pos("dependency-pos"),
	attr_dependency_head("dependency-head"),
	attr_dependency_dependent("dependency-dependent")
    {
      rule_reduce1 = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(1, vocab_type::X)));
      rule_reduce2 = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(2, vocab_type::X)));
    }
    
    typedef utils::bit_vector<1024> coverage_type;
    typedef boost::fusion::tuple<const coverage_type*, int, int> state_type;
    
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<state_type, hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<state_type>,
				    std::allocator<std::pair<state_type, hypergraph_type::id_type> > > state_set_type;
#else
    typedef sgi::hash_map<state_type, hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<state_type>,
			  std::allocator<std::pair<state_type, hypergraph_type::id_type> > > state_set_type;
#endif

#ifdef HAVE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_set<coverage_type, boost::hash<coverage_type>, std::equal_to<coverage_type>,
				    std::allocator<coverage_type > > coverage_set_type;
#else
    typedef sgi::hash_set<coverage_type, boost::hash<coverage_type>, std::equal_to<coverage_type>,
			  std::allocator<coverage_type > > coverage_set_type;
#endif
    
    typedef std::deque<state_type, std::allocator<state_type> > queue_type;
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > terminal_map_type;
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      terminals.clear();
      terminals.resize(lattice.size());
      
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
	  
	  terminals[pos] = node_id;
	}
      }
      
      queue.clear();
      states.clear();
      coverages.clear();
      
      coverage_type __coverage_goal;
      for (size_t i = 0; i != lattice.size(); ++ i)
	__coverage_goal.set(i);
      
      const coverage_type* coverage_start = coverage_vector(coverage_type()).first;
      const coverage_type* coverage_goal  = coverage_vector(__coverage_goal).first;
      
      states[state_type(coverage_start, 0, 0)] = hypergraph_type::invalid;
      queue.push_back(state_type(coverage_start, 0, 0));
      
      hypergraph_type::edge_type::node_set_type tails(2);
      
      // breadth first search
      while (! queue.empty()) {
	const state_type state = queue.front();
	queue.pop_front();
	
	const coverage_type* coverage = boost::fusion::get<0>(state);
	const int head                = boost::fusion::get<1>(state);
	const int range               = boost::fusion::get<2>(state);
	
	const hypergraph_type::id_type node_prev = states[state];
	
	// we need to restrict our range...
	const int first = coverage->select(1, false);
	const int last  = lattice.size();
	for (int i = first; i != last; ++ i) 
	  if (! coverage->test(i)) {
	    coverage_type __coverage_new(*coverage);
	    __coverage_new.set(i);
	    
	    const coverage_type* coverage_new = coverage_vector(__coverage_new).first;
	    
	    if (coverage_new == coverage_goal) {
	      if (! graph.is_valid())
		graph.goal = graph.add_node().id;
	      
	      if (node_prev == hypergraph_type::invalid) {
		tails.back() = terminals[i];
		hypergraph_type::edge_type& edge = graph.add_edge(tails.begin() + 1, tails.end());
		edge.rule = rule_reduce1;
		edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(head);
		edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(i + 1);
		
		graph.connect_edge(edge.id, graph.goal);
	      } else {
		tails.front() = node_prev;
		tails.back() = terminals[i];
		hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		edge.rule = rule_reduce2;
		edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(head);
		edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(i + 1);
		
		graph.connect_edge(edge.id, graph.goal);
	      }
	      
	      // since it is almost the goal, we do not enumerate them!
	      break;
	      
	    } else {
	      // we need to consider two cases... one use the previous-head as our new nead or use dependent as our new head...
	      
	      if (head && i >= range) {
		const state_type state_next(coverage_new, head, i + 1);
		
		std::pair<state_set_type::iterator, bool> result = states.insert(std::make_pair(state_next, 0));
		if (result.second)
		  result.first->second = graph.add_node().id;
		
		if (node_prev == hypergraph_type::invalid) {
		  tails.back() = terminals[i];
		  hypergraph_type::edge_type& edge = graph.add_edge(tails.begin() + 1, tails.end());
		  edge.rule = rule_reduce1;
		  edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(head);
		  edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(i + 1);
		  
		  graph.connect_edge(edge.id, result.first->second);
		} else {
		  tails.front() = node_prev;
		  tails.back() = terminals[i];
		  hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		  edge.rule = rule_reduce2;
		  edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(head);
		  edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(i + 1);
		  
		  graph.connect_edge(edge.id, result.first->second);
		}
		
		if (result.second)
		  queue.push_back(state_next);
	      }
	      
	      const state_type state_next(coverage_new, i + 1, 0);
	      
	      std::pair<state_set_type::iterator, bool> result = states.insert(std::make_pair(state_next, 0));
	      if (result.second)
		result.first->second = graph.add_node().id;
		
	      if (node_prev == hypergraph_type::invalid) {
		tails.back() = terminals[i];
		hypergraph_type::edge_type& edge = graph.add_edge(tails.begin() + 1, tails.end());
		edge.rule = rule_reduce1;
		edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(head);
		edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(i + 1);
		  
		graph.connect_edge(edge.id, result.first->second);
	      } else {
		tails.front() = node_prev;
		tails.back() = terminals[i];
		hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
		edge.rule = rule_reduce2;
		edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(head);
		edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(i + 1);
		  
		graph.connect_edge(edge.id, result.first->second);
	      }
		
	      if (result.second)
		queue.push_back(state_next);
	    }
	  }
      }

      if (graph.is_valid())
	graph.topologically_sort();
    }

    std::pair<const coverage_type*, bool> coverage_vector(const coverage_type& coverage)
    {
      std::pair<coverage_set_type::iterator, bool> result = coverages.insert(coverage);
      
      return std::make_pair(&(*result.first), result.second);
    }
    
  private:
    const attribute_type attr_dependency_pos;
    const attribute_type attr_dependency_head;
    const attribute_type attr_dependency_dependent;

    queue_type        queue;
    state_set_type    states;
    coverage_set_type coverages;
    terminal_map_type terminals;
    
    rule_ptr_type rule_reduce1;
    rule_ptr_type rule_reduce2;
  };

  inline
  void compose_dependency_top_down(const Lattice& lattice, HyperGraph& graph)
  {
    ComposeDependencyTopDown composer;
    composer(lattice, graph);
  }
};

#endif
