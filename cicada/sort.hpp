// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SORT__HPP__
#define __CICADA__SORT__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>

#include <utils/hashmurmur.hpp>

#include <google/dense_hash_set>

namespace cicada
{
  struct TopologicallySort
  {
    // DFS to locale the ordering of hypergraph.
    // The ordering is the same as the post-traversal order of a tree

    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    enum color_type {
      white,
      gray,
      black
    };

    struct dfs_type
    {
      id_type node;
      int edge;
      int tail;
      
      dfs_type(const id_type& _node, const int& _edge, const int& _tail) 
	: node(_node), edge(_edge), tail(_tail) {}
    };
    
    typedef std::vector<int, std::allocator<int> > reloc_set_type;
    typedef std::vector<color_type, std::allocator<color_type> > color_set_type;
    typedef std::vector<dfs_type, std::allocator<dfs_type> > stack_type;

    struct no_filter_edge
    {
      bool operator()(const edge_type& edge) const
      {
	return false;
      }
    };

    struct filter_edge
    {
      std::vector<bool, std::allocator<bool> > removed;

      filter_edge(size_t size) : removed(size, false) {}
      
      bool operator()(const edge_type& edge) const
      {
	return removed[edge.id];
      }
    };
    
    template <typename Filter>
    void operator()(const hypergraph_type& x, hypergraph_type& sorted, Filter filter, const bool validate=true)
    {
      typedef google::dense_hash_set<id_type, utils::hashmurmur<size_t>, std::equal_to<id_type> > id_set_type;
      
      sorted.clear();
      
      if (x.goal == hypergraph_type::invalid)
	return;

      id_set_type edges_cycle;
      edges_cycle.set_empty_key(id_type(-1));

      reloc_set_type reloc_node(x.nodes.size(), -1);
      reloc_set_type reloc_edge(x.edges.size(), -1);
      color_set_type color(x.nodes.size(), white);
      stack_type stack;
      
      stack.reserve(x.nodes.size());
      stack.push_back(dfs_type(x.goal, 0, 0));
      
      int node_count = 0;
      int edge_count = 0;
      
      while (! stack.empty()) {
	const dfs_type& dfs = stack.back();
	id_type node_id = dfs.node;
	size_t pos_edge = dfs.edge;
	size_t pos_tail = dfs.tail;
	
	stack.pop_back();
	
	const node_type* curr_node = &(x.nodes[node_id]);
	
	while (pos_edge != curr_node->edges.size()) {
	  const edge_type& curr_edge = x.edges[curr_node->edges[pos_edge]];
	  
	  if (pos_tail == curr_edge.tails.size() || filter(curr_edge)) {
	    // reach end: proceed to the next edge with pos_tail initialized to the first tail
	    ++ pos_edge;
	    pos_tail = 0;
	    continue;
	  }
	  
	  const id_type tail_node = curr_edge.tails[pos_tail];
	  const color_type tail_color = color[tail_node];
	  
	  switch (tail_color) {
	  case white:
	    ++ pos_tail;
	    stack.push_back(dfs_type(node_id, pos_edge, pos_tail));
	    
	    node_id = tail_node;
	    curr_node = &(x.nodes[node_id]);
	    
	    color[node_id] = gray;
	    pos_edge = 0;
	    pos_tail = 0;
	    
	    break;
	  case black:
	    ++ pos_tail;
	    break;
	  case gray:
	    // cycle detected...
	    // we will force cutting this cycle!
	    ++ pos_tail;

	    edges_cycle.insert(curr_edge.id);
	    
#if 0
	    {
	      std::cerr << "backtrack: " << *curr_edge.rule << std::endl;
	      stack_type::const_reverse_iterator siter_end = stack.rend();
	      for (stack_type::const_reverse_iterator siter = stack.rbegin(); siter != siter_end; ++ siter)
		std::cerr << "backtrack: " << *(x.edges[x.nodes[siter->node].edges[siter->edge]].rule) << std::endl;
	    }
	    
	    throw std::runtime_error("detected cycle!: " + boost::lexical_cast<std::string>(*curr_edge.rule));
#endif
	    
	    break;
	  }
	}
	
	for (size_t i = 0; i != curr_node->edges.size(); ++ i)
	  if (! filter(x.edges[curr_node->edges[i]]))
	    reloc_edge[curr_node->edges[i]] = edge_count ++;
	
	color[node_id] = black;
	reloc_node[node_id] = node_count ++;
      }

#if 0
      if (! edges_cycle.empty()) {
	std::cerr << "cycle detected!" << std::endl;
	std::cerr << x << std::endl;
      }
#endif
      
      // sorted graph!
      sorted.clear();
      
      // construct edges...
      for (size_t i = 0; i != reloc_edge.size(); ++ i)
	if (reloc_edge[i] >= 0) {
	  const edge_type& edge_old = x.edges[i];
	  
	  const id_type edge_id = sorted.edges.size();
	  
	  sorted.edges.push_back(edge_old);
	  
	  edge_type& edge_new = sorted.edges.back();
	  
	  edge_new.id = edge_id;
	  
	  edge_new.head = reloc_node[edge_new.head];
	  edge_type::node_set_type::iterator niter_end = edge_new.tails.end();
	  for (edge_type::node_set_type::iterator niter = edge_new.tails.begin(); niter != niter_end; ++ niter)
	    *niter = reloc_node[*niter];
	  
	  reloc_edge[i] = edge_id;
	}
      
      // construct reverse node-map ...
      reloc_set_type reloc_map_node(node_count, -1);
      for (size_t i = 0; i != x.nodes.size(); ++ i)
	if (reloc_node[i] >= 0)
	  reloc_map_node[reloc_node[i]] = i;

      id_set_type nodes_empty;
      nodes_empty.set_empty_key(id_type(-1));
      
      for (size_t i = 0; i != reloc_map_node.size(); ++ i) {
	const node_type& node_old = x.nodes[reloc_map_node[i]];
	node_type& node_new = sorted.add_node();
	
	node_type::edge_set_type::const_iterator eiter_end = node_old.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node_old.edges.begin(); eiter != eiter_end; ++ eiter)
	  if (reloc_edge[*eiter] >= 0)
	    node_new.edges.push_back(reloc_edge[*eiter]);
	
	if (node_new.edges.empty())
	  nodes_empty.insert(node_new.id);
      }
      
      sorted.goal = sorted.nodes.size() - 1;
      
      if ((! nodes_empty.empty() && validate) || ! edges_cycle.empty()) {
	hypergraph_type sorted_new;
	filter_edge filter(sorted.edges.size());
	
	
	id_set_type::const_iterator eiter_end = edges_cycle.end();
	for (id_set_type::const_iterator eiter = edges_cycle.begin(); eiter != eiter_end; ++ eiter)
	  if (reloc_edge[*eiter] >= 0)
	    filter.removed[reloc_edge[*eiter]] = true;
	
	if (! nodes_empty.empty() && validate)
	  for (typename hypergraph_type::edge_set_type::const_iterator eiter = sorted.edges.begin(); eiter != sorted.edges.end(); ++ eiter) {
	    const edge_type& edge = *eiter;
	    
	    typename edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	    for (typename edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	      if (nodes_empty.find(*titer) != nodes_empty.end()) {
		filter.removed[edge.id] = true;
		break;
	      }
	  }
	
	operator()(sorted, sorted_new, filter, validate);
	
	sorted.swap(sorted_new);
      }
    }
  };
  
  template <typename Filter>
  inline
  void topologically_sort(const HyperGraph& source, HyperGraph& target, Filter filter, const bool validate=true)
  {
    TopologicallySort()(source, target, filter, validate);
  }
  
  inline
  void topologically_sort(const HyperGraph& source, HyperGraph& target, const bool validate=true)
  {
    TopologicallySort()(source, target, TopologicallySort::no_filter_edge(), validate);
  }
  
  inline
  void topologically_sort(HyperGraph& graph)
  {
    HyperGraph x;
    topologically_sort(graph, x);
    graph.swap(x);
  }
};

#endif
