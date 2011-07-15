// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__REMOVE_BOS_EOS__HPP__
#define __CICADA__REMOVE_BOS_EOS__HPP__ 1

#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <utility>

#include <cicada/lattice.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sort.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/mathop.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/bithack.hpp>

namespace cicada
{
  
  struct RemoveBosEos
  {
    // bos_eos removal, but mainly concerns confusion-network lattice...
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Lattice    lattice_type;
    typedef HyperGraph hypergraph_type;
    typedef Vocab      vocab_type;

    typedef hypergraph_type::rule_type rule_type;
    
    typedef lattice_type::symbol_type      symbol_type;
    typedef lattice_type::feature_set_type feature_set_type;
    
    typedef std::pair<int, feature_set_type> bos_eos_type;
    typedef std::set<bos_eos_type, std::less<bos_eos_type>, std::allocator<bos_eos_type> > closure_type;
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
    
    void operator()(const lattice_type& source, lattice_type& target)
    {
      
      target.clear();
      
      if (source.empty()) return;
      
      // first, compute closure
      
      closure_set_type closure(source.size() + 1);
      
      for (int state = source.size() - 1; state >= 0; -- state) {
	lattice_type::arc_set_type::const_iterator aiter_end = source[state].end();
	for (lattice_type::arc_set_type::const_iterator aiter = source[state].begin(); aiter != aiter_end; ++ aiter) 
	  if (aiter->label == vocab_type::BOS || aiter->label == vocab_type::EOS) {
	    const int last = state + aiter->distance;
	    
	    closure[state].insert(bos_eos_type(last, aiter->features));
	    
	    closure_type::const_iterator citer_end = closure[last].end();
	    for (closure_type::const_iterator citer = closure[last].begin(); citer != citer_end; ++ citer)
	      closure[state].insert(bos_eos_type(citer->first, citer->second + aiter->features));
	  }
      }

      // then, compute graph with closure...
      
      graph_type removed(source.size() + 1);
      backptr_set_type backptr(source.size() + 1);
      
      for (size_t state = 0; state != source.size(); ++ state) {
	lattice_type::arc_set_type::const_iterator aiter_end = source[state].end();
	for (lattice_type::arc_set_type::const_iterator aiter = source[state].begin(); aiter != aiter_end; ++ aiter) 
	  if (aiter->label != vocab_type::BOS && aiter->label != vocab_type::EOS) {
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
	      if (niter->label != vocab_type::BOS && niter->label != vocab_type::EOS) {
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
      
      // sort removed wrt edge...
      removed.sort_edge();
      
      // pruning ...
      position_set_type positions_removed(removed.nodes.size(), -1);
      const size_type  num_nodes_removed = dfs(removed, 0, positions_removed);

      if (positions_removed.back() != 0)
	std::cerr << "WARNING: DFS resulted in wrong lattice (1st step)" << std::endl;
      
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

      if (positions_transposed.back() != 0)
	std::cerr << "WARNING: DFS resulted in wrong lattice (2nd step)" << std::endl;
      
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

    typedef utils::simple_vector<int, std::allocator<int> > index_set_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > rhs_set_type;
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tail_set_type;
    
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > bos_eos_set_type;
    typedef std::vector<bos_eos_set_type, std::allocator<bos_eos_set_type> > bos_eos_map_type;

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
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      
      target = source;

      if (! target.is_valid()) return;

      rhs_set_type rhs;
      tail_set_type tails;
      
      bos_eos_map_type bos_eoss(target.nodes.size());
      removed_type     removed(target.edges.size(), false);
      size_type        bos_eos_remove = 0;
      size_type        bos_eos_new = 0;
      
      hypergraph_type::node_set_type::iterator niter_end = target.nodes.end();
      for (hypergraph_type::node_set_type::iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	const hypergraph_type::node_type& node_source = source.nodes[node.id];
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  hypergraph_type::edge_type& edge = target.edges[*eiter];
	  
	  const rule_type& rule = *edge.rule;
	  
	  if (rule.rhs.size() == 1 && (rule.rhs.front() == vocab_type::BOS || rule.rhs.front() == vocab_type::EOS)) {
	    // we will mark this as deleted,
	    // and keep bit-vector indicating the node.id has a edge with bos_eos...
	    
	    bos_eoss[node.id].push_back(edge.id);
	    removed[edge.id] = true;
	    ++ bos_eos_remove;
	  } else {
	    rhs.clear();
	    rule_type::symbol_set_type::const_iterator siter_end = rule.rhs.end();
	    for (rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter) {
	      if (*siter != vocab_type::BOS && *siter != vocab_type::EOS)
		rhs.push_back(*siter);
	    }
	    
	    if (rhs.size() != rule.rhs)
	      edge.rule = rule_type::create(rule_type(rule.lhs, rhs.begin(), rhs.end()));
	    
	    index_set_type j_ends(edge.tails.size(), 0);
	    index_set_type j(edge.tails.size(), 0);
	    
	    bool found_bos_eos = false;
	    
	    for (size_type i = 0; i != edge.tails.size(); ++ i) {
	      found_bos_eos |= ! bos_eoss[edge.tails[i]].empty();
	      j_ends[i] = utils::bithack::branch(bos_eoss[edge.tails[i]].empty(), size_type(0), bos_eoss[edge.tails[i]].size() + 1);
	    }
	    
	    if (! found_bos_eos) continue;
	    
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
		  
		  if (j[antecedent_index] > 0 && j_ends[antecedent_index] > 0) {
		    const hypergraph_type::edge_type& edge_antecedent = target.edges[bos_eoss[edge.tails[antecedent_index]][j[antecedent_index] - 1]];
		    
		    features += edge_antecedent.features;
		    // edge_antcedent is an bos_eos rule without tails!
		  } else {
		    tails.push_back(edge.tails[antecedent_index]);
		    rhs.push_back(riter->non_terminal());
		  }
		} else
		  rhs.push_back(*riter);
	      }

	      // this new bos_eos-rule should be deleted again!
	      if (tails.empty() && rhs.empty()) {
		rhs.push_back(vocab_type::BOS);
		++ bos_eos_new;
	      }
	      
	      hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	      edge_new.rule = rule_type::create(rule_type(edge.rule->lhs, rhs.begin(), rhs.end()));
	      edge_new.features = features;
	      edge_new.attributes = edge.attributes;
	      
	      target.connect_edge(edge_new.id, edge.head);
	      
	      // proceed to the next...
	      size_type index = 0;
	      for (/**/; index != j.size(); ++ index) 
		if (j_ends[index] != 0){
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
      if (bos_eos_remove) {
	removed.resize(target.edges.size(), false);
	
	hypergraph_type sorted;
	topologically_sort(target, sorted, filter_edge(removed), true);
	target.swap(sorted);
      }
      
      // further removal...
      if (target.is_valid() && bos_eos_new) {
	hypergraph_type removed;
	operator()(target, removed);
	target.swap(removed);
      }
    }
  };
  
  inline
  void remove_bos_eos(Lattice& lattice)
  {
    RemoveBosEos __remover;
    Lattice __lattice;
    __remover(lattice, __lattice);
    lattice.swap(__lattice);
  }
  
  inline
  void remove_bos_eos(const Lattice& lattice, Lattice& removed)
  {
    RemoveBosEos __remover;
    __remover(lattice, removed);
  }

  inline
  void remove_bos_eos(HyperGraph& graph)
  {
    RemoveBosEos __remover;
    HyperGraph removed;
    __remover(graph, removed);
    graph.swap(removed);
  }
  
  inline
  void remove_bos_eos(const HyperGraph& graph, HyperGraph& removed)
  {
    RemoveBosEos __remover;
    __remover(graph, removed);
  }
  
};


#endif
