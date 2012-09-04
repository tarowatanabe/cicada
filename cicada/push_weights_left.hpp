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
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef std::vector<id_type, std::allocator<id_type> > edge_list_type;
    typedef std::vector<edge_list_type, std::allocator<edge_list_type> > edge_set_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;

      edge_set_type outgoing(target.nodes.size());
      edge_set_type leaning(target.nodes.size());
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = target.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = target.edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = *eiter;
	
	if (edge.tails.empty()) continue;
	
	int tail_pos = 0;
	rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) 
	  if (riter->is_non_terminal()) {
	    const int __non_terminal_index = riter->non_terminal_index();
	    const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, tail_pos, __non_terminal_index - 1);
	    
	    if (antecedent_index == 0)
	      leaning[edge.tails[tail_pos]].push_back(edge.id);
	    
	    outgoing[edge.tails[tail_pos]].push_back(edge.id);
	    
	    ++ tail_pos;
	  }
      }
      
      hypergraph_type::node_set_type::const_reverse_iterator niter_end = target.nodes.rend();
      for (hypergraph_type::node_set_type::const_reverse_iterator niter = target.nodes.rbegin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	// this should not happen, though..
	if (node.edges.empty()) continue;
	if (outgoing[node.id].empty()) continue;
	
	const edge_list_type& accumulate = (leaning[node.id].empty() ? outgoing[node.id] : leaning[node.id]);
	
	if (accumulate.size() == 1) {
	  
	  if (! target.edges[accumulate.front()].features.empty()) {
	    const feature_set_type intersected(target.edges[accumulate.front()].features);
	    
	    node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	      target.edges[*eiter].features += intersected;
	    
	    edge_list_type::const_iterator oiter_end = outgoing[node.id].end();
	    for (edge_list_type::const_iterator oiter = outgoing[node.id].begin(); oiter != oiter_end; ++ oiter)
	      target.edges[*oiter].features -= intersected;
	  }
	} else {
	  feature_set_type intersected(target.edges[accumulate.front()].features);
	    
	  edge_list_type::const_iterator oiter_end = accumulate.end();
	  for (edge_list_type::const_iterator oiter = accumulate.begin() + 1; oiter != oiter_end; ++ oiter)
	    intersected.intersect(target.edges[*oiter].features);

	  if (! intersected.empty()) {
	    edge_list_type::const_iterator oiter_end = outgoing[node.id].end();
	    for (edge_list_type::const_iterator oiter = outgoing[node.id].begin(); oiter != oiter_end; ++ oiter)
	      target.edges[*oiter].features -= intersected;
	      
	    node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	      target.edges[*eiter].features += intersected;
	  }
	}
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
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef std::vector<id_type, std::allocator<id_type> > edge_list_type;
    typedef std::vector<edge_list_type, std::allocator<edge_list_type> > edge_set_type;

    typedef std::vector<bool, std::allocator<bool> > root_set_type;
    
    typedef hypergraph_type::edge_type::node_set_type tails_type;
    typedef std::vector<tails_type, std::allocator<tails_type> > tails_set_type;
    
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
    
    typedef std::vector<int, std::allocator<int> > position_set_type;
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
      tails_set_type tails(target.edges.size() + 1);
      edge_set_type incoming(target.nodes.size() + 1);
      edge_set_type outgoing(target.nodes.size());
      edge_set_type neutral(target.nodes.size());
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = target.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = target.edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = *eiter;

	if (edge.tails.empty()) {
	  neutral[edge.head].push_back(edge.id);
	  continue;
	}
	
	outgoing[edge.head].push_back(edge.id);

	roots[edge.head] = false;
	tails[edge.id] = edge.tails;
	
	int tail_pos = 0;
	rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) 
	  if (riter->is_non_terminal()) {
	    const int __non_terminal_index = riter->non_terminal_index();
	    const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, tail_pos, __non_terminal_index - 1);
	    
	    if (antecedent_index == 0) {
	      incoming[edge.tails[tail_pos]].push_back(edge.id);
	      tails[edge.id][tail_pos] = edge.head;
	    } else
	      outgoing[edge.tails[tail_pos]].push_back(edge.id);
	    
	    ++ tail_pos;
	  }
      }
      
      // add a pseudo root with a hyperedge...
      {
	std::vector<id_type, std::allocator<id_type> > tails_root;
        
	for (id_type node = 0; node != roots.size(); ++ node)
	  if (roots[node])
	    tails_root.push_back(node);
	
	tails.back() = tails_type(tails_root.begin(), tails_root.end());
	incoming.back().push_back(tails.size() - 1);
      }
      
      // topologically sort...
      position_set_type positions(incoming.size(), -1);
      color_set_type    color(incoming.size(), white);
      stack_type        stack;
      
      stack.reserve(incoming.size());
      stack.push_back(dfs_type(incoming.size() - 1, 0, 0));
      
      int node_count = 0;
      
      while (! stack.empty()) {
	const dfs_type& dfs = stack.back();
	id_type node_id = dfs.node;
	size_t pos_edge = dfs.edge;
	size_t pos_tail = dfs.tail;
	
	stack.pop_back();
       
	const edge_list_type* curr_node = &incoming[node_id];
       
	while (pos_edge != curr_node->size()) {
	  const tails_type& curr_edge = tails[curr_node->operator[](pos_edge)];
         
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
	    curr_node = &incoming[node_id];
	    
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
	positions[node_count ++] = node_id;
      }
      
      position_set_type::const_iterator piter_end = positions.end() - 1;
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	if (*piter < 0) continue;
	
	const hypergraph_type::node_type& node = target.nodes[*piter];
	

#if 0	
	// we will collect weights from incoming and neutral, then, distribute to outging...
	
	if (incoming[node.id].empty() && neutral[node.id].empty()) continue;
	if (outgoing[node.id].empty()) continue;
	
	feature_set_type intersected;
	
	if (! incoming[node.id].empty()) {
	  intersected = target.edges[incoming[node.id].front()].features;
	  
	  edge_list_type::const_iterator iiter_end = incoming[node.id].end();
	  for (edge_list_type::const_iterator iiter = incoming[node.id].begin() + 1; iiter != iiter_end; ++ iiter)
	    intersected.intersect(target.edges[*iiter].features);
	}
	
	if (! neutral[node.id].empty()) {
	  if (incoming[node.id].empty())
	    intersected = target.edges[neutral[node.id].front()].features;
	  else
	    intersected.intersect(target.edges[neutral[node.id].front()].features);
	  
	  edge_list_type::const_iterator niter_end = neutral[node.id].end();
	  for (edge_list_type::const_iterator niter = neutral[node.id].begin() + 1; niter != niter_end; ++ niter)
	    intersected.intersect(target.edges[*niter].features);
	}
	
	if (! intersected.empty()) {
	  edge_list_type::const_iterator iiter_end = incoming[node.id].end();
	  for (edge_list_type::const_iterator iiter = incoming[node.id].begin(); iiter != iiter_end; ++ iiter)
	    target.edges[*iiter].features -= intersected;
	  
	  edge_list_type::const_iterator niter_end = neutral[node.id].end();
	  for (edge_list_type::const_iterator niter = neutral[node.id].begin(); niter != niter_end; ++ niter)
	    target.edges[*niter].features -= intersected;
	  
	  edge_list_type::const_iterator oiter_end = outgoing[node.id].end();
	  for (edge_list_type::const_iterator oiter = outgoing[node.id].begin(); oiter != oiter_end; ++ oiter)
	    target.edges[*oiter].features += intersected;
	}
#endif
	
#if 1
	// TODO!
	// we will collect weights from outging, then distribute to incoming...
	// if no outgoing, we will treat neutral as outging. otherwise neutral is incoming
	
	const edge_list_type& accumulate = (incoming[node.id].empty() ? neutral[node.id] : incoming[node.id]);
	
	if (accumulate.empty()) continue;
	
	feature_set_type intersected;
	
	if (! outgoing[node.id].empty()) {
	  intersected = target.edges[outgoing[node.id].front()].features;
	  
	  edge_list_type::const_iterator oiter_end = outgoing[node.id].end();
	  for (edge_list_type::const_iterator oiter = outgoing[node.id].begin() + 1; oiter != oiter_end; ++ oiter)
	    intersected.intersect(target.edges[*oiter].features);
	}
	
	if (! incoming[node.id].empty() && ! neutral[node.id].empty()) {
	  if (outgoing[node.id].empty())
	    intersected = target.edges[neutral[node.id].front()].features;
	  else
	    intersected.intersect(target.edges[neutral[node.id].front()].features);
	  
	  edge_list_type::const_iterator niter_end = neutral[node.id].end();
	  for (edge_list_type::const_iterator niter = neutral[node.id].begin() + 1; niter != niter_end; ++ niter)
	    intersected.intersect(target.edges[*niter].features);
	}
	
	if (! intersected.empty()) {
#if 0
	  std::cerr << "propagated: " << node.id
		    << " incoming: " << incoming[node.id].size()
		    << " outgoing: " << outgoing[node.id].size()
		    << " neutral: " << neutral[node.id].size()
		    << std::endl;
	  std::cerr << intersected;
#endif

	  edge_list_type::const_iterator oiter_end = outgoing[node.id].end();
	  for (edge_list_type::const_iterator oiter = outgoing[node.id].begin(); oiter != oiter_end; ++ oiter)
	    target.edges[*oiter].features -= intersected;
	  
	  if (! incoming[node.id].empty()) {
	    edge_list_type::const_iterator niter_end = neutral[node.id].end();
	    for (edge_list_type::const_iterator niter = neutral[node.id].begin(); niter != niter_end; ++ niter)
	      target.edges[*niter].features -= intersected;
	  }
	  
	  edge_list_type::const_iterator aiter_end = accumulate.end();
	  for (edge_list_type::const_iterator aiter = accumulate.begin(); aiter != aiter_end; ++ aiter)
	    target.edges[*aiter].features += intersected;
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
