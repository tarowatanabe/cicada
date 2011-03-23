// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__REMOVE_UNKNOWN__HPP__
#define __CICADA__REMOVE_UNKNOWN__HPP__ 1

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/sort.hpp>

#include <utils/bithack.hpp>
#include <utils/simple_vector.hpp>

namespace cicada
{
  struct RemoveUnknown
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    typedef rule_type::symbol_type     symbol_type;
    typedef rule_type::symbol_set_type symbol_set_type;

    typedef std::vector<bool, std::allocator<bool> > unknown_type;
    typedef std::vector<bool, std::allocator<bool> > removed_type;

    typedef utils::simple_vector<int, std::allocator<int> > index_set_type;
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tail_set_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > rhs_type;
    
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
      if (! source.is_valid()) {
	target.clear();
	return;
      }

      target = source;

      removed_type removed(target.edges.size(), false);
      unknown_type unknown(target.edges.size(), false);

      // first, check whether it is unk...
      hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	// this should not happen, though..
	if (node.edges.empty()) continue;
	
	const hypergraph_type::edge_type& edge = target.edges[node.edges.front()];
	const symbol_type& lhs = edge.rule->lhs;
	
	if (! is_unknown(lhs)) continue;
	
	unknown[node.id] = true;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	  removed[*eiter] = true;
      }
      
      tail_set_type tails;
      rhs_type      rhs;
      
      for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const node_type& node = *niter;
	
	const size_type edges_size = node.edges.size();
	for (size_type e = 0; e != edges_size; ++ e) {
	  const hypergraph_type::edge_type& edge = target.edges[node.edges[e]];
	  
	  // search for antecedent nodes, and seek the unknown label..
	  // if found, try merge! 
	  // it is like apply-exact to form new edges....
	  
	  index_set_type j_ends(edge.tails.size(), 0);
	  index_set_type j(edge.tails.size(), 0);
	  
	  bool found_unknown = false;
	  
	  for (size_type i = 0; i != edge.tails.size(); ++ i) {
	    found_unknown |= unknown[edge.tails[i]];
	    j_ends[i] = utils::bithack::branch(unknown[edge.tails[i]], target.nodes[edge.tails[i]].edges.size(), size_type(0));
	  }
	  
	  if (! found_unknown) continue;
	  
	  removed[edge.id] = true;
	  
	  // TODO: we do not care index in categories...
	  
	  for (;;) {
	    tails.clear();
	    rhs.clear();
	    
	    feature_set_type features = edge.features;
	    
	    bool invalid = false;
	    size_type i = 0;
	    rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	    for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
	      if (riter->is_non_terminal()) {
		if (j_ends[i] > 0) {
		  const hypergraph_type::node_type& node_antecedent = target.nodes[edge.tails[i]];
		  const hypergraph_type::edge_type& edge_antecedent = target.edges[node_antecedent.edges[j[i]]];
		  
		  features += edge_antecedent.features;
		  
		  hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge_antecedent.tails.end();
		  for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge_antecedent.tails.begin(); titer != titer_end; ++ titer) {
		    invalid |= unknown[*titer];
		    tails.push_back(*titer);
		  }
		  
		  rule_type::symbol_set_type::const_iterator aiter_end = edge_antecedent.rule->rhs.end();
		  for (rule_type::symbol_set_type::const_iterator aiter = edge_antecedent.rule->rhs.begin(); aiter != aiter_end; ++ aiter)
		    rhs.push_back(aiter->is_non_terminal() ? aiter->non_terminal() : *aiter);
		} else {
		  tails.push_back(edge.tails[i]);
		  rhs.push_back(riter->non_terminal());
		}
		
		++ i;
	      } else
		rhs.push_back(*riter);
	    
	    if (! invalid) {
	      hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	      edge_new.rule = rule_type::create(rule_type(edge.rule->lhs, rhs.begin(), rhs.end()));
	      edge_new.features = features;
	      edge_new.attributes = edge.attributes;
	      
	      target.connect_edge(edge_new.id, edge.head);

	      if (removed.size() < target.edges.size())
		removed.resize(target.edges.size(), false);
	      
	      removed[edge_new.id] = unknown[edge.head];
	    }
	    
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
      
      removed.resize(target.edges.size(), false);
      
      //target.topologically_sort();
      
      hypergraph_type sorted;
      topologically_sort(target, sorted, filter_edge(removed), true);
      target.swap(sorted);
    }

    bool is_unknown(const symbol_type& symbol)
    {
      if (! symbol.is_non_terminal()) return;
      
      const utils::piece piece = symbol.non_terminal_strip();
      
      return (piece.size() >= 7 && piece.substr(0, 7) == "UNKNOWN");
    }
    
  };
  
  inline
  void remove_unknown(HyperGraph& graph)
  {
    RemoveUnknown remover;
    
    HyperGraph graph_removed;
    remover(graph, graph_removed);
    graph.swap(graph_removed);
  }
  
  inline
  void remove_unknown(const HyperGraph& source, HyperGraph& target)
  {
    RemoveUnknown remover;
    
    remover(source, target);
  }

};

#endif
