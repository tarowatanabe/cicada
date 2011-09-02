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

#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/bithack.hpp>
#include <utils/indexed_trie.hpp>

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
    
    typedef boost::fusion::tuple<int, int, int> span_type; // parent, span
    typedef utils::indexed_trie<span_type, utils::hashmurmur<size_t>, std::equal_to<span_type>, std::allocator<span_type> > stack_type;
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_set_type;
    typedef std::deque<stack_type::id_type, std::allocator<stack_type::id_type> > queue_type;
    typedef std::vector<bool, std::allocator<bool> > queued_type;
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      terminals.clear();
      
      // initialize actives by axioms... (terminals)

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
      
      stack.clear();
      nodes.clear();
      queue.clear();
      
      // initialize stack..
      const stack_type::id_type id = stack.push(stack.root(), span_type(0, 1, lattice.size() + 1));
      if (id != 0)
	throw std::runtime_error("we assume id is zero!");
      nodes.resize(1, hypergraph_type::invalid);
      queued.resize(1, true);
      queue.push_back(id);
      
      hypergraph_type::edge_type::node_set_type tails(2);
      
      // breadth first search
      while (! queue.empty()) {
	const stack_type::id_type state = queue.front();
	queue.pop_front();
	
	const int head  = boost::fusion::get<0>(stack[state]);
	const int first = boost::fusion::get<1>(stack[state]);
	const int last  = boost::fusion::get<2>(stack[state]);

	const stack_type::id_type state_prev = stack.pop(state);
	const hypergraph_type::id_type node_prev = nodes[state];
	
	for (int dependent = first; dependent != last; ++ dependent) {
	  stack_type::id_type state_next = state_prev;
	  
	  if (dependent + 1 != last) {
	    state_next = stack.push(state_next, span_type(dependent, dependent + 1, last));
	    if (state_next >= nodes.size())
	      nodes.resize(state_next + 1, hypergraph_type::invalid);
	    if (nodes[state_next] == hypergraph_type::invalid)
	      nodes[state_next] = graph.add_node().id;
	  }
	  
	  if (first != dependent) {
	    state_next = stack.push(state_next, span_type(dependent, first, dependent));
	    if (state_next >= nodes.size())
	      nodes.resize(state_next + 1, hypergraph_type::invalid);
	    if (nodes[state_next] == hypergraph_type::invalid)
	      nodes[state_next] = graph.add_node().id;
	  }
	  
	  hypergraph_type::id_type node_parent;
	  
	  // we reached goal
	  if (state_next == stack.root()) {
	    if (graph.goal == hypergraph_type::invalid)
	      graph.goal = graph.add_node().id;
	    
	    node_parent = graph.goal;
	  } else
	    node_parent = nodes[state_next];
	  
	  if (node_prev == hypergraph_type::invalid) {
	    tails.back() = terminals[dependent - 1];
	    hypergraph_type::edge_type& edge = graph.add_edge(tails.begin() + 1, tails.end());
	    edge.rule = rule_reduce1;
	    edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(head);
	    edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(dependent);
	    
	    graph.connect_edge(edge.id, node_parent);
	  } else {
	    tails.front() = node_prev;
	    tails.back()  = terminals[dependent - 1];
	    hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
	    edge.rule = rule_reduce2;
	    edge.attributes[attr_dependency_head]      = attribute_set_type::int_type(head);
	    edge.attributes[attr_dependency_dependent] = attribute_set_type::int_type(dependent);
	    
	    graph.connect_edge(edge.id, node_parent);
	  }

	  if (state_next != stack.root()) {
	    
	    // push state_next into queue... if already queues, igore!
	    if (state_next >= queued.size())
	      queued.resize(state_next + 1, false);
	    
	    if (! queued[state_next]) {
	      queue.push_back(state_next);
	      queued[state_next] = true;
	    }
	  }
	}
      }

      if (graph.is_valid())
	graph.topologically_sort();
    }
    
  private:
    const attribute_type attr_dependency_pos;
    const attribute_type attr_dependency_head;
    const attribute_type attr_dependency_dependent;

    queue_type    queue;
    queued_type   queued;
    stack_type    stack;
    node_set_type nodes;
    node_set_type terminals;
    
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
