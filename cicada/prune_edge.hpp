// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PRUNE_EDGE__HPP__
#define __CICADA__PRUNE_EDGE__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sort.hpp>

namespace cicada
{
  template <typename Function>
  struct PruneEdge
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    typedef Function function_type;
    
    typedef typename function_type::value_type weight_type;
    
    typedef std::vector<bool, std::allocator<bool> > removed_type;

    struct filter_pruned
    {
      const removed_type& removed;
      
      filter_pruned(const removed_type& __removed) : removed(__removed) {}
      
      template <typename Edge>
      bool operator()(const Edge& edge) const
      {
	return removed[edge.id];
      }
    };


    PruneEdge(const function_type& __function,
	      const size_t __size,
	      const bool __validate=true)
      : function(__function),
	size(__size),
	validate(__validate) {}
    
    typedef std::vector<weight_type, std::allocator<weight_type> > inside_type;
    typedef std::vector<weight_type, std::allocator<weight_type> > posterior_type;

    typedef std::pair<weight_type, id_type> value_type;
    typedef std::vector<value_type, std::allocator<value_type> > sorted_type;
    
    template <typename Tp>
    struct greater_first
    {
      bool operator()(const Tp& x, const Tp& y) const { return x.first > y.first; }
    };
  
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      
      target.clear();
      
      if (! source.is_valid())
	return;
      
      inside_type    inside(source.nodes.size());
      posterior_type posterior(source.edges.size());
      
      inside_outside(source, inside, posterior, function, function);
      
      removed_type removed(source.edges.size(), false);
      sorted_type  sorted;
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const node_type& node = *niter;
	
	if (node.edges.size() <= size) continue;
	
	sorted.clear();
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter)
	  sorted.push_back(std::make_pair(posterior[*eiter], *eiter));
	
	std::nth_element(sorted.begin(), sorted.begin() + size, sorted.end(), greater_first<value_type>());
	
	const weight_type threshold = sorted[size].first;
	
	typename sorted_type::const_iterator siter = sorted.begin();
	typename sorted_type::const_iterator siter_end = sorted.end();
	typename sorted_type::const_iterator siter_last = siter + size;
	
	bool found_equal = false;
	for (/**/; siter != siter_last; ++ siter)
	  found_equal |= (siter->first == threshold);
	
	if (found_equal) {
	  for (/**/; siter != siter_end; ++ siter)
	    removed[siter->second] = (siter->first != threshold);
	} else {
	  for (/**/; siter != siter_end; ++ siter)
	    removed[siter->second] = true;
	}
      }
      
      topologically_sort(source, target, filter_pruned(removed), validate);
    }

    const function_type& function;
    const size_t size;
    const bool validate;
  };
  
  
  template <typename Function>
  inline
  void prune_edge(const HyperGraph& source, HyperGraph& target, const Function& func, const size_t size, const bool validate=true)
  {
    PruneEdge<Function> __prune(func, size, validate);
    
    __prune(source, target);
  }
  
  template <typename Function>
  inline
  void prune_edge(HyperGraph& source, const Function& func, const size_t size, const bool validate=true)
  {
    PruneEdge<Function> __prune(func, size, validate);

    HyperGraph target;
    
    __prune(source, target);
    
    source.swap(target);
  }

};

#endif
