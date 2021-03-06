// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__REMOVE_SYMBOL__HPP__
#define __CICADA__REMOVE_SYMBOL__HPP__ 1

#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <utility>

#include <cicada/lattice.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sort_topologically.hpp>
#include <cicada/feature_vector_linear.hpp>

#include <utils/mathop.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/bithack.hpp>
#include <utils/small_vector.hpp>

namespace cicada
{
  
  struct RemoveSymbol
  {
    // epsilon removal, but mainly concerns confusion-network lattice...
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Lattice    lattice_type;
    typedef HyperGraph hypergraph_type;
    typedef Vocab      vocab_type;

    typedef hypergraph_type::rule_type rule_type;
    
    typedef lattice_type::symbol_type      symbol_type;
    typedef lattice_type::feature_set_type feature_set_type;

    typedef FeatureVectorLinear<double> feature_linear_type;
    
    typedef std::pair<int, feature_linear_type> epsilon_type;

    typedef std::set<epsilon_type, std::less<epsilon_type>, std::allocator<epsilon_type> > closure_type;
    typedef std::vector<closure_type, std::allocator<closure_type> > closure_set_type;

    typedef std::vector<int, std::allocator<int> > position_set_type;
    typedef std::vector<int, std::allocator<int> > stack_type;

    typedef std::set<int, std::less<int>, std::allocator<int> > backptr_type;
    typedef std::vector<backptr_type, std::allocator<backptr_type> > backptr_set_type;

    enum color_type {
      white,
      gray,
      black
    };

    typedef std::vector<color_type, std::allocator<color_type> >  color_set_type;
    
    template <typename __Graph>
    size_type dfs(const __Graph& graph, const int start, position_set_type& positions)
    {
      size_type pos = 0;
      
      stack_type stack;
      color_set_type color(positions.size(), white);
      
      stack.push_back(start);
      
      while (! stack.empty()) {
	const int state = stack.back();
	switch (color[state]) {
	case white:
	  color[state] = gray;
	  {
	    typename __Graph::const_iterator aiter_end = graph[state].end();
	    for (typename __Graph::const_iterator aiter = graph[state].begin(); aiter != aiter_end; ++ aiter) {
	      const int state_next = graph.antecedent(state, *aiter);
	      if (color[state_next] == white)
		stack.push_back(state_next);
	      // otherwise, cycle detected...
	    }
	  }
	  break;
	case gray:
	  color[state] = black;
	  positions[state] = pos ++;

	  stack.pop_back();
	  break;
	case black:
	  stack.pop_back();
	  break;
	}
      }
      
      return pos;
    }
    
    struct Graph
    {
      struct edge_type
      {
	int id;
	symbol_type label;
	feature_set_type features;
	int tail;

	edge_type(const symbol_type& __label, const feature_set_type& __features, const int __tail)
	  : id(-1), label(__label), features(__features), tail(__tail) {}
      };
      
      struct node_type
      {
	typedef std::vector<int, std::allocator<int> > edge_set_type;
	
	node_type() : edges() {}
	
	edge_set_type edges;
      };

      typedef utils::chunk_vector<edge_type, 4096 / sizeof(edge_type), std::allocator<edge_type> > edge_set_type;
      typedef std::vector<node_type, std::allocator<node_type> > node_set_type;
      
      typedef node_type::edge_set_type::const_iterator const_iterator;
      
      const node_type::edge_set_type& operator[](size_type pos) const
      {
	return nodes[pos].edges;
      }
      
      size_type antecedent(const size_type pos, const int id) const
      {
	return edges[id].tail;
      }
      
      bool valid(size_type pos) const
      {
	return pos == nodes.size() - 1 || ! nodes[pos].edges.empty();
      }

      struct less_edge
      {
	less_edge(const Graph& __graph) : graph(__graph) {}
	
	const Graph& graph;
	
	bool operator()(const int& x, const int& y) const
	{
	  const edge_type& edge1 = graph.edges[x];
	  const edge_type& edge2 = graph.edges[y];
	  
	  return edge1.tail < edge2.tail;
#if 0
	  return (edge1.tail < edge2.tail
		  || (edge1.tail == edge2.tail
		      && (edge1.label < edge2.label
			  || (edge1.label == edge2.label
			      && edge1.features < edge2.features))));
#endif
	}
      };
      
      void sort_edge()
      {
	node_set_type::iterator niter_end = nodes.end();
	for (node_set_type::iterator niter = nodes.begin(); niter != niter_end; ++ niter)
	  std::sort(niter->edges.begin(), niter->edges.end(), less_edge(*this));

#if 0
	std::set<int, less_edge, std::allocator<int> > sorted(less_edge(*this));

	node_set_type::iterator niter_end = nodes.end();
	for (node_set_type::iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	  sorted.clear();
	  sorted.insert(niter->edges.begin(), niter->edges.end());
	  
	  niter->edges.assign(sorted.begin(), sorted.end());
	}
#endif
      }
      
      edge_set_type edges;
      node_set_type nodes;

      Graph() : edges(), nodes() {}
      Graph(size_type size) : edges(), nodes(size) {}

      void clear() { edges.clear(); nodes.clear(); }

      edge_type& add_edge(const symbol_type& label, const feature_set_type& features, const int tail) 
      {
	const int id = edges.size();
	
	edges.push_back(edge_type(label, features, tail));
	edges.back().id = id;
	
	return edges.back();
      }
      
    };
    

    struct LatticeGraph
    {
      typedef lattice_type::arc_set_type::const_iterator const_iterator;
      
      LatticeGraph(const lattice_type& __lattice) : lattice(__lattice), arcs() {}
      
      const lattice_type::arc_set_type& operator[](size_type pos) const
      {
	return (pos == lattice.size() ? arcs : lattice[pos]);
      }
      
      size_type antecedent(const size_type pos, const lattice_type::arc_type& arc) const
      {
	return pos + arc.distance;
      }
      
      bool valid(size_type pos) const
      {
	return pos == lattice.size() || ! lattice[pos].empty();
      }
      
      const lattice_type& lattice;
      const lattice_type::arc_set_type arcs;
    };

    
    typedef Graph        graph_type;
    typedef LatticeGraph lattice_graph_type;
    
    template <typename RemoveSymbol>
    void operator()(const lattice_type& source, lattice_type& target, RemoveSymbol remove_symbol)
    {
      
      target.clear();
      
      if (source.empty()) return;
      
      // first, compute closure
      
      closure_set_type closure(source.size() + 1);
      
      for (int state = source.size() - 1; state >= 0; -- state) {
	lattice_type::arc_set_type::const_iterator aiter_end = source[state].end();
	for (lattice_type::arc_set_type::const_iterator aiter = source[state].begin(); aiter != aiter_end; ++ aiter) 
	  if (remove_symbol(aiter->label)) {
	    const int last = state + aiter->distance;
	    
	    closure[state].insert(epsilon_type(last, feature_linear_type(aiter->features)));
	    
	    closure_type::const_iterator citer_end = closure[last].end();
	    for (closure_type::const_iterator citer = closure[last].begin(); citer != citer_end; ++ citer)
	      closure[state].insert(epsilon_type(citer->first, feature_linear_type(feature_set_type(citer->second) + aiter->features)));
	  }
      }

      // then, compute graph with closure...
      
      graph_type removed(source.size() + 1);
      backptr_set_type backptr(source.size() + 1);
      
      for (size_t state = 0; state != source.size(); ++ state) {
	lattice_type::arc_set_type::const_iterator aiter_end = source[state].end();
	for (lattice_type::arc_set_type::const_iterator aiter = source[state].begin(); aiter != aiter_end; ++ aiter) 
	  if (! remove_symbol(aiter->label)) {
	    const int state_next = state + aiter->distance;
	    
	    graph_type::edge_type& edge = removed.add_edge(aiter->label, aiter->features, state_next);
	    removed.nodes[state].edges.push_back(edge.id);
	    backptr[state_next].insert(state);
	  }
	
	closure_type::const_iterator citer_end = closure[state].end();
	for (closure_type::const_iterator citer = closure[state].begin(); citer != citer_end; ++ citer)
	  if (citer->first != static_cast<int>(source.size())) {
	    lattice_type::arc_set_type::const_iterator niter_end = source[citer->first].end();
	    for (lattice_type::arc_set_type::const_iterator niter = source[citer->first].begin(); niter != niter_end; ++ niter) 
	      if (! remove_symbol(niter->label)) {
		const int state_next = citer->first + niter->distance;

		graph_type::edge_type& edge = removed.add_edge(niter->label, feature_set_type(citer->second) + niter->features, state_next);
		removed.nodes[state].edges.push_back(edge.id);
		backptr[state_next].insert(state);
	      }
	  } else {
	    // collect edges to the end, indicated by closure...
	    
	    backptr_type::const_iterator biter_end = backptr[state].end();
	    for (backptr_type::const_iterator biter = backptr[state].begin(); biter != biter_end; ++ biter) {
	      const int state_prev = *biter;
	      
	      const size_type arc_size = removed.nodes[state_prev].edges.size();
	      for (size_type i = 0; i != arc_size; ++ i) {
		const graph_type::edge_type& arc = removed.edges[removed.nodes[state_prev].edges[i]];
		
		if (arc.tail == static_cast<int>(state)) {
		  const int state_next = source.size();
		  
		  graph_type::edge_type& edge = removed.add_edge(arc.label, arc.features + feature_set_type(citer->second), state_next);
		  removed.nodes[state_prev].edges.push_back(edge.id);
		  backptr[state_next].insert(state_prev);
		}
	      }
	    }
	  }
      }
      
      // sort removed wrt edge...
      removed.sort_edge();
      
      // pruning ...
      position_set_type positions_removed(removed.nodes.size(), -1);
      const size_type  num_nodes_removed = dfs(removed, 0, positions_removed);
      
      const int goal_removed = positions_removed.back();

      if (goal_removed < 0)
	throw std::runtime_error("DFS resulted in wrong lattice (1st step)");
      
      // after dfs, positons_removed is numberd by post-traversal order... thus,
      // we can automatically transpose the graph!

      graph_type transposed(num_nodes_removed - goal_removed);
      
      for (size_t i = 0; i != positions_removed.size(); ++ i)
	if (positions_removed[i] >= goal_removed) {
	  graph_type::node_type::edge_set_type::const_iterator eiter_end = removed.nodes[i].edges.end();
	  for (graph_type::node_type::edge_set_type::const_iterator eiter = removed.nodes[i].edges.begin(); eiter != eiter_end; ++ eiter) {
	    const graph_type::edge_type& edge = removed.edges[*eiter];
	    
	    if (positions_removed[edge.tail] >= goal_removed) {
	      graph_type::edge_type& edge_new = transposed.add_edge(edge.label, edge.features, positions_removed[i] - goal_removed);
	      
	      transposed.nodes[positions_removed[edge.tail] - goal_removed].edges.push_back(edge_new.id);
	    }
	  }
	}
      
      position_set_type positions_transposed(transposed.nodes.size(), -1);
      const size_type num_nodes_transposed = dfs(transposed, 0, positions_transposed);

      const int goal_transposed = positions_transposed.back();

      if (goal_transposed < 0)
	throw std::runtime_error("DFS resulted in wrong lattice (2nd step)");
      
      // after dfs, positons_transposed is numberd by post-traversal order... thus,
      // we can automatically transpose the graph... combined with the previous transposition,
      // we can uncover original pruned graph!
      
      target.resize(num_nodes_transposed - 1 - goal_transposed);
      
      for (size_t i = 0; i != positions_transposed.size(); ++ i)
	if (positions_transposed[i] >= goal_transposed) {
	  graph_type::node_type::edge_set_type::const_iterator eiter_end = transposed.nodes[i].edges.end();
	  for (graph_type::node_type::edge_set_type::const_iterator eiter = transposed.nodes[i].edges.begin(); eiter != eiter_end; ++ eiter) {
	    const graph_type::edge_type& edge = transposed.edges[*eiter];
	    
	    if (positions_transposed[edge.tail] >= goal_transposed) {
	      const int first = positions_transposed[edge.tail] - goal_transposed;
	      const int last  = positions_transposed[i] - goal_transposed;
	      
	      target[first].push_back(lattice_type::arc_type(edge.label, edge.features, last - first));
	    }
	  }
	}
      
      target.initialize_distance();
    }

    typedef utils::small_vector<int, std::allocator<int> > index_set_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > rhs_set_type;
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tail_set_type;
    
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > epsilon_set_type;
    typedef std::vector<epsilon_set_type, std::allocator<epsilon_set_type> > epsilon_map_type;

    typedef std::vector<bool, std::allocator<bool> > removed_type;
    
    struct filter_edge
    {
      filter_edge(const removed_type& __removed) : removed(__removed) {}
      
      bool operator()(const hypergraph_type::edge_type& edge) const
      {
        return removed[edge.id];
      }
      
      const removed_type& removed;
     };
    
    template <typename RemoveSymbol>
    void operator()(const hypergraph_type& source, hypergraph_type& target, RemoveSymbol remove_symbol)
    {
      
      target = source;

      if (! target.is_valid()) return;

      rhs_set_type rhs;
      tail_set_type tails;
      
      epsilon_map_type epsilons(target.nodes.size());
      removed_type     removed(target.edges.size(), false);
      size_type        epsilon_remove = 0;
      
      hypergraph_type::node_set_type::iterator niter_end = target.nodes.end();
      for (hypergraph_type::node_set_type::iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	const size_t edge_size = node.edges.size();
	for (size_t e = 0; e != edge_size; ++ e) {
	  hypergraph_type::edge_type& edge = target.edges[node.edges[e]];
	  
	  const rule_type& rule = *edge.rule;
	  
	  if (edge.tails.empty()) {
	    // we will check terminals...
	    if (rule.rhs.size() == 1) {
	      if (remove_symbol(rule.rhs.front())) {
		epsilons[node.id].push_back(edge.id);
		removed[edge.id] = true;
		++ epsilon_remove;
	      }
	    } else {
	      rhs.clear();
	      rule_type::symbol_set_type::const_iterator siter_end = rule.rhs.end();
	      for (rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter)
		if (! remove_symbol(*siter))
		  rhs.push_back(*siter);
	      
	      if (rhs.size() != rule.rhs.size()) {
		if (rhs.empty()) {
		  epsilons[node.id].push_back(edge.id);
		  removed[edge.id] = true;
		  ++ epsilon_remove;
		  
		  rhs.push_back(vocab_type::EPSILON);
		}
		
		edge.rule = rule_type::create(rule_type(rule.lhs, rhs.begin(), rhs.end()));
	      }
	    }
	  } else {
	    rhs.clear();
	    rule_type::symbol_set_type::const_iterator siter_end = rule.rhs.end();
	    for (rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter)
	      if (! remove_symbol(*siter))
		rhs.push_back(*siter);
	    
	    if (rhs.size() != rule.rhs.size())
	      edge.rule = rule_type::create(rule_type(rule.lhs, rhs.begin(), rhs.end()));
	      
	    index_set_type j_ends(edge.tails.size(), 0);
	    index_set_type j(edge.tails.size(), 0);
	    
	    bool found_epsilon = false;
	    
	    for (size_type i = 0; i != edge.tails.size(); ++ i) {
	      found_epsilon |= ! epsilons[edge.tails[i]].empty();
	      j_ends[i] = epsilons[edge.tails[i]].size();
	    }
	    
	    if (! found_epsilon) continue;
	    
	    for (;;) {
	      tails.clear();
	      rhs.clear();
	      
	      feature_set_type features = edge.features;
	      
	      int non_terminal_pos = 0;
	      rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	      for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) {
		if (riter->is_non_terminal()) {
		  const int __non_terminal_index = riter->non_terminal_index();
		  const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
		  ++ non_terminal_pos;
		  
		  if (j_ends[antecedent_index]) {
		    const hypergraph_type::edge_type& edge_antecedent = target.edges[epsilons[edge.tails[antecedent_index]][j[antecedent_index]]];
		    
		    features += edge_antecedent.features;
		    // edge_antcedent is an epsilon rule without tails!
		  } else {
		    tails.push_back(edge.tails[antecedent_index]);
		    rhs.push_back(riter->non_terminal());
		  }
		} else
		  rhs.push_back(*riter);
	      }
	      
	      // this new epsilon-rule should be deleted...
	      bool epsilon_new = false;
	      if (tails.empty() && rhs.empty()) {
		rhs.push_back(vocab_type::EPSILON);
		epsilon_new = true;
	      }
	      
	      hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	      edge_new.rule = rule_type::create(rule_type(edge.rule->lhs, rhs.begin(), rhs.end()));
	      edge_new.features = features;
	      edge_new.attributes = edge.attributes;
	      
	      target.connect_edge(edge_new.id, edge.head);

	      if (epsilon_new) {
		epsilons[node.id].push_back(edge_new.id);
		
		if (edge_new.id >= removed.size())
		  removed.resize(edge_new.id + 1, false);
		removed[edge_new.id] = true;
		++ epsilon_remove;
	      }
	      
	      // proceed to the next...
	      size_type index = 0;
	      for (/**/; index != j.size(); ++ index) 
		if (j_ends[index]) {
		  ++ j[index];
		  if (j[index] < j_ends[index]) break;
		  j[index] = 0;
		}
	      
	      // finished!
	      if (index == j.size()) break;
	    }
	  }
	}
      }
      
      // topologically sort...
      if (epsilon_remove) {
	removed.resize(target.edges.size(), false);
	
	hypergraph_type sorted;
	topologically_sort(target, sorted, filter_edge(removed), true);
	target.swap(sorted);
      }
    }
  };
  
  template <typename RemoveSymbol>
  inline
  void remove_symbol(Lattice& lattice, RemoveSymbol rs)
  {
    RemoveSymbol __remover;
    Lattice __lattice;
    __remover(lattice, __lattice, rs);
    lattice.swap(__lattice);
  }
  
  template <typename RemoveSymbol>
  inline
  void remove_symbol(const Lattice& lattice, Lattice& removed, RemoveSymbol rs)
  {
    RemoveSymbol __remover;
    __remover(lattice, removed, rs);
  }
  
  template <typename RemoveSymbol>
  inline
  void remove_symbol(HyperGraph& graph, RemoveSymbol rs)
  {
    RemoveSymbol __remover;
    HyperGraph removed;
    __remover(graph, removed, rs);
    graph.swap(removed);
  }
  
  template <typename RemoveSymbol>
  inline
  void remove_symbol(const HyperGraph& graph, HyperGraph& removed, RemoveSymbol rs)
  {
    RemoveSymbol __remover;
    __remover(graph, removed, rs);
  }
  
};


#endif
