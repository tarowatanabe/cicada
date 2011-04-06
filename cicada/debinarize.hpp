// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__DEBINARIZE__HPP__
#define __CICADA__DEBINARIZE__HPP__ 1

#include <stdexcept>
#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/sort.hpp>

#include <utils/bithack.hpp>
#include <utils/simple_vector.hpp>

namespace cicada
{
  struct Debinarize
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::symbol_type      symbol_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef std::vector<bool, std::allocator<bool> > binarized_type;
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
      // debinarization by stripping off the ^ from syntactic categories...
      
      // bottom-up topological order to find binarised antecedents...
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;
      
      removed_type removed(target.edges.size(), false);
      binarized_type binarized(target.nodes.size(), false);
      
      // first, check whether it is binarized!
      hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	// this should not happen, though..
	if (node.edges.empty()) continue;
	
	const hypergraph_type::edge_type& edge = target.edges[node.edges.front()];
	const symbol_type& lhs = edge.rule->lhs;
	
	if (lhs.non_terminal_strip().find('^') != symbol_type::piece_type::npos()) {
	  binarized[node.id] = true;
	  
	  hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	    removed[*eiter] = true;
	}
      }

      tail_set_type tails;
      rhs_type      rhs;
      
      for (hypergraph_type::node_set_type::const_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;

	const size_type edges_size = node.edges.size();
	for (size_type e = 0; e != edges_size; ++ e) {
	  const hypergraph_type::edge_type& edge = target.edges[node.edges[e]];
	  
	  // search for antecedent nodes, and seek the binarized label..
	  // if found, try merge! 
	  // it is like apply-exact to form new edges....
	  
	  //std::cerr << "rule: " << *edge.rule << std::endl;
	  
	  index_set_type j_ends(edge.tails.size(), 0);
	  index_set_type j(edge.tails.size(), 0);
	  
	  bool found_binarized = false;
	  
	  for (size_type i = 0; i != edge.tails.size(); ++ i) {
	    found_binarized |= binarized[edge.tails[i]];
	    j_ends[i] = utils::bithack::branch(binarized[edge.tails[i]], target.nodes[edge.tails[i]].edges.size(), size_type(0));
	  }
	  
	  if (! found_binarized) continue;
	  
	  removed[edge.id] = true;
	  
	  // TODO: we do not care index in categories...
	  // This is potentially a hard work in that we need to adjust index...
	  
	  for (;;) {
	    tails.clear();
	    rhs.clear();
	    
	    feature_set_type features = edge.features;
	    
	    bool invalid = false;
	    int non_terminal_pos = 0;
	    rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	    for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
	      if (riter->is_non_terminal()) {
		const int __non_terminal_index = riter->non_terminal_index();
		const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
		++ non_terminal_pos;
		
		if (j_ends[antecedent_index] > 0) {
		  const hypergraph_type::node_type& node_antecedent = target.nodes[edge.tails[antecedent_index]];
		  const hypergraph_type::edge_type& edge_antecedent = target.edges[node_antecedent.edges[j[antecedent_index]]];
		  
		  features += edge_antecedent.features;
		  
		  // special care is reqiured for gran-antecedents by converting indices....
		  
		  hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge_antecedent.tails.end();
		  for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge_antecedent.tails.begin(); titer != titer_end; ++ titer) {
		    invalid |= binarized[*titer];
		    tails.push_back(*titer);
		  }
		  
		  rule_type::symbol_set_type::const_iterator aiter_end = edge_antecedent.rule->rhs.end();
		  for (rule_type::symbol_set_type::const_iterator aiter = edge_antecedent.rule->rhs.begin(); aiter != aiter_end; ++ aiter)
		    rhs.push_back(aiter->is_non_terminal() ? aiter->non_terminal() : *aiter);
		  
		} else {
		  tails.push_back(edge.tails[antecedent_index]);
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
	      
	      removed[edge_new.id] = binarized[edge.head];
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
    
  };


  
  inline
  void debinarize(const HyperGraph& source, HyperGraph& target)
  {
    Debinarize debinarizer;
    
    debinarizer(source, target);
  }
  
  inline
  void debinarize(HyperGraph& source)
  {
    HyperGraph target;
    debinarize(source, target);
    source.swap(target);
  }

};

#endif
