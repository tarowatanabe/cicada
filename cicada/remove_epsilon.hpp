// -*- mode: c++ -*-

#ifndef __CICADA__REMOVE_EPSILON__HPP__
#define __CICADA__REMOVE_EPSILON__HPP__ 1

#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <utility>

#include <cicada/lattice.hpp>
#include <cicada/vocab.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/mathop.hpp>
#include <utils/chunk_vector.hpp>

namespace cicada
{
  
  struct RemoveEpsilon
  {
    // epsilon removal, but mainly concerns confusion-network lattice...
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Lattice lattice_type;
    typedef Vocab   vocab_type;
    
    typedef lattice_type::symbol_type      symbol_type;
    typedef lattice_type::feature_set_type feature_set_type;
    
    typedef std::pair<int, feature_set_type> epsilon_type;
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
      
      const lattice_type& lattice;
      const lattice_type::arc_set_type arcs;
    };

    
    typedef Graph        graph_type;
    typedef LatticeGraph lattice_graph_type;
    
    void operator()(const lattice_type& source, lattice_type& target)
    {
      
      target.clear();
      
      if (source.empty()) return;
      
      // first, compute closure
      
      closure_set_type closure(source.size() + 1);
      
      for (int state = source.size() - 1; state >= 0; -- state) {
	lattice_type::arc_set_type::const_iterator aiter_end = source[state].end();
	for (lattice_type::arc_set_type::const_iterator aiter = source[state].begin(); aiter != aiter_end; ++ aiter) 
	  if (aiter->label == vocab_type::EPSILON) {
	    const int last = state + aiter->distance;
	    
	    closure[state].insert(epsilon_type(last, aiter->features));
	    
	    closure_type::const_iterator citer_end = closure[last].end();
	    for (closure_type::const_iterator citer = closure[last].begin(); citer != citer_end; ++ citer)
	      closure[state].insert(epsilon_type(citer->first, citer->second + aiter->features));
	  }
      }

      // then, compute graph with closure...
      
      graph_type removed(source.size() + 1);
      backptr_set_type backptr(source.size() + 1);
      
      for (size_t state = 0; state != source.size(); ++ state) {
	lattice_type::arc_set_type::const_iterator aiter_end = source[state].end();
	for (lattice_type::arc_set_type::const_iterator aiter = source[state].begin(); aiter != aiter_end; ++ aiter) 
	  if (aiter->label != vocab_type::EPSILON) {
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
	      if (niter->label != vocab_type::EPSILON) {
		const int state_next = citer->first + niter->distance;

		graph_type::edge_type& edge = removed.add_edge(niter->label, citer->second + niter->features, state_next);
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
		  
		  graph_type::edge_type& edge = removed.add_edge(arc.label, arc.features + citer->second, state_next);
		  removed.nodes[state_prev].edges.push_back(edge.id);
		  backptr[state_next].insert(state_prev);
		}
	      }
	    }
	  }
      }
      
      // pruning ...
      position_set_type positions_removed(removed.nodes.size(), -1);
      const size_type  num_nodes_removed = dfs(removed, 0, positions_removed);
      
      // after dfs, positons_removed is numberd by post-traversal order... thus,
      // we can automatically transpose the graph!

      graph_type transposed(num_nodes_removed);
      
      for (size_t i = 0; i != positions_removed.size(); ++ i)
	if (positions_removed[i] >= 0) {
	  graph_type::node_type::edge_set_type::const_iterator eiter_end = removed.nodes[i].edges.end();
	  for (graph_type::node_type::edge_set_type::const_iterator eiter = removed.nodes[i].edges.begin(); eiter != eiter_end; ++ eiter) {
	    const graph_type::edge_type& edge = removed.edges[*eiter];
	    
	    if (positions_removed[edge.tail] >= 0) {
	      graph_type::edge_type& edge_new = transposed.add_edge(edge.label, edge.features, positions_removed[i]);
	      
	      transposed.nodes[positions_removed[edge.tail]].edges.push_back(edge_new.id);
	    }
	  }
	}
      
      position_set_type positions_transposed(transposed.nodes.size(), -1);
      const size_type num_nodes_transposed = dfs(transposed, 0, positions_transposed);
      
      // after dfs, positons_transposed is numberd by post-traversal order... thus,
      // we can automatically transpose the graph... combined with the previous transposition,
      // we can uncover original pruned graph!
      
      target.resize(num_nodes_transposed - 1);
      
      for (size_t i = 0; i != positions_transposed.size(); ++ i)
	if (positions_transposed[i] >= 0) {
	  graph_type::node_type::edge_set_type::const_iterator eiter_end = transposed.nodes[i].edges.end();
	  for (graph_type::node_type::edge_set_type::const_iterator eiter = transposed.nodes[i].edges.begin(); eiter != eiter_end; ++ eiter) {
	    const graph_type::edge_type& edge = transposed.edges[*eiter];
	    
	    if (positions_transposed[edge.tail] >= 0) {
	      const int first = positions_transposed[edge.tail];
	      const int last  = positions_transposed[i];
	      
	      target[first].push_back(lattice_type::arc_type(edge.label, edge.features, last - first));
	    }
	  }
	}
      
      
      target.initialize_distance();
    }
  };
  
  inline
  void remove_epsilon(Lattice& lattice)
  {
    RemoveEpsilon __remover;
    Lattice __lattice;
    __remover(lattice, __lattice);
    lattice.swap(__lattice);
  }
  

  inline
  void remove_epsilon(const Lattice& lattice, Lattice& removed)
  {
    RemoveEpsilon __remover;
    __remover(lattice, removed);
  }
  
};


#endif
