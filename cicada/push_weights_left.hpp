// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// one of the major problem is the interaction from other tails...
// HOW TO SOLVE THIS???
// use of inside/outside???
//

#ifndef __CICADA__PUSH_WEIGHTS_LEFT__HPP__
#define __CICADA__PUSH_WEIGHTS_LEFT__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>

#include <utils/bithack.hpp>

namespace cicada
{
  
  struct PushWeightsLeft
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type   id_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    typedef std::vector<feature_set_type, std::allocator<feature_set_type> > feature_map_type;
    
    typedef std::vector<id_type, std::allocator<id_type> > edge_list_type;
    typedef std::vector<edge_list_type, std::allocator<edge_list_type> > edge_set_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      
    }
  };

  
#if 0
  struct PushWeightsLeft
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type   id_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    typedef std::vector<feature_set_type, std::allocator<feature_set_type> > feature_map_type;
    
    typedef std::vector<id_type, std::allocator<id_type> > edge_list_type;
    typedef std::vector<edge_list_type, std::allocator<edge_list_type> > edge_set_type;

    typedef std::vector<id_type, std::allocator<id_type> > tail_set_type;

    typedef std::vector<bool, std::allocator<bool> > root_set_type;
    
    enum color_type {
      white,
      gray,
      black
    };
    
    typedef std::vector<color_type, std::allocator<color_type> >  color_set_type;
    
    typedef std::vector<int, std::allocator<int> > position_set_type;
    typedef std::vector<int, std::allocator<int> > stack_type;

    typedef std::vector<bool, std::allocator<bool> > processed_type;

    //
    // we will construct a lattice structure, then, propagate weights...
    // for source-node, we will distribute weights assigned to each node...
    //
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;
      
      edge_set_type incoming(target.nodes.size() + 1);
      tail_set_type tails(target.edges.size());
      root_set_type roots(target.nodes.size(), true);
      
      // visit in a topological order, then, compute a list of incoming edges
      {
	hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
	for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	  const hypergraph_type::node_type& node = *niter;
	  
	  hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	    const hypergraph_type::edge_type& edge = target.edges[*eiter];
	    
	    if (edge.tails.empty()) continue;
	    
	    int tail_pos = 0;
	    rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	    for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) 
	      if (riter->is_non_terminal()) {
		const int __non_terminal_index = riter->non_terminal_index();
		const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, tail_pos, __non_terminal_index - 1);
		
		if (antecedent_index == 0) break;
		
		++ tail_pos;
	      }
	    
	    incoming[edge.tails[tail_pos]].push_back(edge.id);
	    tails[edge.id] = edge.head;
	    roots[edge.head] = false;
	  }
	}
      }
      
      for (id_type node_id = 0; node_id != roots.size(); ++ node_id)
	if (roots[node_id]) {
	  incoming.back().push_back(tails.size());
	  tails.push_back(node_id);
	}

      std::cerr << "incomings: " << incoming.size() << " tails: " << tails.size() << std::endl;
      
      // depth-first-search...
      
      size_type node_pos = 0;
      position_set_type positions(incoming.size(), -1);
      
      color_set_type color(incoming.size(), white);
      stack_type stack;
      
      stack.reserve(incoming.size());
      stack.push_back(incoming.size() - 1);
      
      while (! stack.empty()) {
	const id_type node = stack.back();
	
	switch (color[node]) {
	case white:
	  color[node] = gray;
	  {
	    edge_list_type::const_iterator aiter_end = incoming[node].end();
	    for (edge_list_type::const_iterator aiter = incoming[node].begin(); aiter != aiter_end; ++ aiter) {
	      const id_type next = tails[*aiter];
	      
	      if (color[next] == white)
		stack.push_back(next);
	      // otherwise, cycle detected...
	    }
	  }
	  break;
	case gray:
	  color[node] = black;
	  //positions[node] = node_pos ++;
	  positions[node_pos ++] = node;
	  
	  stack.pop_back();
	  break;
	case black:
	  stack.pop_back();
	  break;
	}
      }
      
      feature_map_type features(target.nodes.size());
      //processed_type processed(positions.size(), false);
      
      position_set_type::const_iterator piter_end = positions.end() - 1;
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	//processed[*piter] = true;

	if (incoming[*piter].empty()) continue;
	
	edge_list_type::const_iterator iiter_end = incoming[*piter].end();
	for (edge_list_type::const_iterator iiter = incoming[*piter].begin(); iiter != iiter_end; ++ iiter) {
	  target.edges[*iiter].features += features[tails[*iiter]];
	  
	  //if (! processed[tails[*iiter]])
	  //  std::cerr << "WARNING: not processed?: " << *piter << " tail: " << tails[*iiter] << std::endl;
	}

	if (incoming[*piter].size() == 1) {
	  features[*piter] = target.edges[incoming[*piter].front()].features;
	  
	  target.edges[incoming[*piter].front()].features.clear();
	} else {
	  feature_set_type intersected(target.edges[incoming[*piter].front()].features);
	  
	  edge_list_type::const_iterator iiter_end = incoming[*piter].end();
	  for (edge_list_type::const_iterator iiter = incoming[*piter].begin() + 1; iiter != iiter_end; ++ iiter)
	    intersected.intersect(target.edges[*iiter].features);
	  
	  if (! intersected.empty()) {
	    features[*piter] = intersected;
	    
	    edge_list_type::const_iterator iiter_end = incoming[*piter].end();
	    for (edge_list_type::const_iterator iiter = incoming[*piter].begin(); iiter != iiter_end; ++ iiter)
	      target.edges[*iiter].features -= intersected;
	  }
	}
      }
      
      hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	   hypergraph_type::edge_type& edge = target.edges[*eiter];
	   
	   if (edge.tails.empty())
	     edge.features += features[node.id];
	}
      }
    }
  };
#endif

#if 0
  struct PushWeightsLeft
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type   id_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    typedef std::vector<feature_set_type, std::allocator<feature_set_type> > feature_map_type;
    
    typedef std::vector<id_type, std::allocator<id_type> > edge_list_type;
    typedef std::vector<edge_list_type, std::allocator<edge_list_type> > edge_set_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;
      
      edge_set_type incoming(target.nodes.size());
      
      // visit in a topological order, then, compute a list of incoming edges
      {
	hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
	for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	  const hypergraph_type::node_type& node = *niter;
	  
	  hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	    const hypergraph_type::edge_type& edge = target.edges[*eiter];
	    
	    if (edge.tails.empty()) continue;
	    
	    int tail_pos = 0;
	    rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	    for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) 
	      if (riter->is_non_terminal()) {
		const int __non_terminal_index = riter->non_terminal_index();
		const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, tail_pos, __non_terminal_index - 1);
		
		if (antecedent_index == 0) break;
		
		++ tail_pos;
	      }
	    
	    incoming[edge.tails[tail_pos]].push_back(edge.id);
	  }
	}
      }
      
      // visit in an inverse of topilogical order, then, propagate weights
      
      hypergraph_type::node_set_type::const_reverse_iterator niter_end = target.nodes.rend();
      for (hypergraph_type::node_set_type::const_reverse_iterator niter = target.nodes.rbegin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	if (incoming[node.id].empty()) continue;

	feature_set_type intersected(target.edges[incoming[node.id].front()].features);
	
	if (incoming[node.id].size() == 1)
	  target.edges[incoming[node.id].front()].features.clear();
	else {
	  edge_list_type::const_iterator iiter_end = incoming[node.id].end();
	  for (edge_list_type::const_iterator iiter = incoming[node.id].begin() + 1; iiter != iiter_end; ++ iiter)
	    intersected.intersect(target.edges[*iiter].features);
	  
	  if (! intersected.empty()) {
	    edge_list_type::const_iterator iiter_end = incoming[node.id].end();
	    for (edge_list_type::const_iterator iiter = incoming[node.id].begin(); iiter != iiter_end; ++ iiter)
	      target.edges[*iiter].features -= intersected;
	  }
	}
	
	// we will also propagate to antecedents
	
	if (! intersected.empty()) {
	  hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	    target.edges[*eiter].features += intersected;
	}
      }
    }
  };
#endif

#if 0
  struct PushWeightsLeft
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type   id_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    typedef std::vector<feature_set_type, std::allocator<feature_set_type> > feature_map_type;
    typedef std::vector<bool, std::allocator<bool> > processed_type;

    typedef std::vector<bool, std::allocator<bool> > root_set_type;

    typedef hypergraph_type::edge_type::node_set_type tail_set_type;
    
    typedef std::vector<tail_set_type, std::allocator<tail_set_type> > edge_set_type;
    
    typedef std::vector<id_type, std::allocator<id_type> > edge_list_type;
    typedef std::vector<edge_list_type, std::allocator<edge_list_type> > node_set_type;
    
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

    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;
      
      root_set_type roots(target.nodes.size(), true);
      node_set_type nodes(target.nodes.size() + 1);
      edge_set_type edges(target.edges.size() + 1);
      
      // visit in a topological order, then, compute a "left-corner graph"
      {
	hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
	for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	  const hypergraph_type::node_type& node = *niter;
	  
	  hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	    const hypergraph_type::edge_type& edge = target.edges[*eiter];
	    
	    if (edge.tails.empty()) {
	      nodes[edge.head].push_back(edge.id);
	      continue;
	    }
	    
	    int tail_pos = 0;
	    rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	    for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) 
	      if (riter->is_non_terminal()) {
		const int __non_terminal_index = riter->non_terminal_index();
		const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, tail_pos, __non_terminal_index - 1);
		
		if (antecedent_index == 0) break;
		
		++ tail_pos;
	      }
	    
	    nodes[edge.tails[tail_pos]].push_back(edge.id);
	    
	    roots[edge.head] = false;
	    
	    edges[edge.id] = edge.tails;
	    edges[edge.id][tail_pos] = edge.head;
	  }
	}
      }
      
      // add a pseudo root with a hyperedge...
      {
	std::vector<id_type, std::allocator<id_type> > tails;
	
	for (id_type node = 0; node != roots.size(); ++ node)
	  if (roots[node])
	    tails.push_back(node);
	
	edges.back() = tail_set_type(tails.begin(), tails.end());
	nodes.back().push_back(edges.size() - 1);
      }
      
      // topologically sort...
      reloc_set_type reloc_node(nodes.size(), -1);
      color_set_type color(nodes.size(), white);
      stack_type stack;
      
      stack.reserve(nodes.size());
      stack.push_back(dfs_type(nodes.size() - 1, 0, 0));
      
      int node_count = 0;
      
      while (! stack.empty()) {
	const dfs_type& dfs = stack.back();
	id_type node_id = dfs.node;
	size_t pos_edge = dfs.edge;
	size_t pos_tail = dfs.tail;
	
	stack.pop_back();
	
	const edge_list_type* curr_node = &nodes[node_id];
	
	while (pos_edge != curr_node->size()) {
	  const tail_set_type& curr_edge = edges[curr_node->operator[](pos_edge)];
	  
	  if (pos_tail == curr_edge.size()) {
	    // reach end: proceed to the next edge with pos_tail initialized to the first tail
	    ++ pos_edge;
	    pos_tail = 0;
	    continue;
	  }
	  
	  const id_type    tail_node  = curr_edge[pos_tail];
	  const color_type tail_color = color[tail_node];
	  
	  switch (tail_color) {
	  case white:
	    ++ pos_tail;
	    stack.push_back(dfs_type(node_id, pos_edge, pos_tail));
	    
	    node_id = tail_node;
	    curr_node = &nodes[node_id];
	    
	    color[node_id] = gray;
	    pos_edge = 0;
	    pos_tail = 0;
	    
	    break;
	  case black:
	    ++ pos_tail;
	    
	    break;
	  case gray:
	    ++ pos_tail;
	    
	    break;
	  }
	}
	
	color[node_id] = black;
	reloc_node[node_count ++] = node_id;
      }
      
      feature_map_type features(target.nodes.size());
      
      reloc_set_type::const_iterator riter_end = reloc_node.end() - 1;
      for (reloc_set_type::const_iterator riter = reloc_node.begin(); riter != riter_end; ++ riter) {
	
	// this should not happen...
	if (*riter < 0) continue;
	if (nodes[*riter].empty()) continue;
	
	// first, propagate features toward edges
	edge_list_type::const_iterator eiter_end = nodes[*riter].end();
	for (edge_list_type::const_iterator eiter = nodes[*riter].begin(); eiter != eiter_end; ++ eiter) {
	  hypergraph_type::edge_type& edge = target.edges[*eiter];
	  
	  tail_set_type::const_iterator titer_end = edges[*eiter].end();
	  for (tail_set_type::const_iterator titer = edges[*eiter].begin(); titer != titer_end; ++ titer)
	    edge.features += features[*titer];
	}
	
	if (roots[*riter]) continue;
	
	// second, compute intersection, and propagate intersected weights.
	
	if (nodes[*riter].size() == 1) {
	  features[*riter] += target.edges[nodes[*riter].front()].features;
	  
	  target.edges[nodes[*riter].front()].features.clear();
	} else {
	  feature_set_type intersected(target.edges[nodes[*riter].front()].features);
	  
	  edge_list_type::const_iterator eiter_end = nodes[*riter].end();
	  for (edge_list_type::const_iterator eiter = nodes[*riter].begin() + 1; eiter != eiter_end; ++ eiter)
	    intersected.intersect(target.edges[*eiter].features);
	  
	  if (! intersected.empty()) {
	    features[*riter] += intersected;
	    
	    edge_list_type::const_iterator eiter_end = nodes[*riter].end();
	    for (edge_list_type::const_iterator eiter = nodes[*riter].begin(); eiter != eiter_end; ++ eiter)
	      target.edges[*eiter].features -= intersected;
	  }
	}
	
#if 0
	if (! features[*riter].empty()) {
	  const hypergraph_type::node_type& node = target.nodes[*riter];
	  
	  hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	    hypergraph_type::edge_type& edge = target.edges[*eiter];
	    
	    if (edge.tails.empty())
	      edge.features += features[*riter];
	  }
	}
#endif
      }
    }
  };
#endif
  
  inline
  void push_weights_left(const HyperGraph& source, HyperGraph& target)
  {
    PushWeightsLeft()(source, target);
  }
  
  inline
  void push_weights_left(HyperGraph& graph)
  {
    HyperGraph x;
    push_weights_left(graph, x);
    graph.swap(x);
  }
};

#endif
