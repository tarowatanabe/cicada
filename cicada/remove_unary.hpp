// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__REMOVE_UNARY__HPP__
#define __CICADA__REMOVE_UNARY__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sort.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/mathop.hpp>

namespace cicada
{
  
  struct RemoveUnary
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef HyperGraph hypergraph_type;
    typedef Vocab      vocab_type;

    typedef hypergraph_type::rule_type rule_type;

    typedef utils::simple_vector<int, std::allocator<int> > index_set_type;
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
      
      removed_type removed(target.edges.size(), false);
      size_type    unary_remove = 0;
      
      hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	const size_t edge_size = node.edges.size();
	for (const size_t e = 0; e != edge_size; ++ e) {
	  const hypergraph_type::edge_type& edge = target.edges[node.edges[e]];
	  
	  const rule_type& rule = *edge.rule;
	  
	  if (edge.tails.size() != 1 || rule.rhs.size() != 1) continue;
	  
	  removed[edge.id] = true;
	  ++ unary_remove;
	  
	  const hypergraph_type::node_type& node_tail = target.nodes[edge.tails.front()];
	  
	  hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_tail.edges.end();
	  for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_tail.edges.begin(); eiter != eiter_end; ++ eiter) {
	    const hypergraph_type::edge_type& edge_antecedent = target.edges[*eiter];
	    const rule_type& rule_antecedent = *edge_antecedent.rule;
	    
	    // ignore unary rules...
	    if (edge_antecedent.tails.size() == 1 && rule.rhs.size() == 1) continue;
	    
	    // create new edge!
	    
	    hypergraph_type::edge_type& edge_new = target.add_edge(edge_antecedent.tails.begin(), edge_antecedent.tails.end());
	    edge_new.rule = rule_type::create(rule_type(rule.lhs, rule_antecedent.rhs));
	    edge_new.attributes = edge.attributes;
	    edge_new.features = edge.features + edge_antecedent.features;
	    
	    target.connect_edge(edge_new.id, edge.head);
	  }
	}
      }
      
      if (unary_remove) {
	removed.resize(target.edges.size(), false);
	
	hypergraph_type sorted;
	topologically_sort(target, sorted, filter_edge(removed), true);
	target.swap(sorted);
      }
    }
  };

  inline
  void remove_unary(HyperGraph& graph)
  {
    RemoveUnary __remover;
    HyperGraph removed;
    __remover(graph, removed);
    graph.swap(removed);
  }
  
  inline
  void remove_unary(const HyperGraph& graph, HyperGraph& removed)
  {
    RemoveUnary __remover;
    __remover(graph, removed);
  }
  
};

#endif
