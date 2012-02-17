// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PRUNE_DENSITY__HPP__
#define __CICADA__PRUNE_DENSITY__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sort_topologically.hpp>
#include <cicada/viterbi.hpp>
#include <cicada/inside_outside.hpp>

namespace cicada
{
  
  template <typename Function>
  struct PruneDensity
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

    struct length_traversal
    {
      typedef int value_type;
      
      template <typename Edge, typename Iterator>
      void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
      {
	yield = 1;
	for (/**/; first != last; ++ first)
	  yield += *first;
      }
    };


    PruneDensity(const function_type& __function,
		 const double __threshold,
		 const bool __validate=true)
      : function(__function),
	threshold(__threshold),
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
      weight_type viterbi_weight;
      int         viterbi_length;
      viterbi(source, viterbi_length, viterbi_weight, length_traversal(), function);
      
      const size_t prune_size = static_cast<size_t>(threshold * viterbi_length);
      
      if (source.edges.size() <= prune_size) {
	target = source;
	return;
      }
      
      inside_type    inside(source.nodes.size());
      posterior_type posterior(source.edges.size());
      
      inside_outside(source, inside, posterior, function, function);
      
      sorted_type sorted(source.edges.size());
      
      for (id_type id = 0; id != source.edges.size(); ++ id)
	sorted[id] = std::make_pair(posterior[id], id);
      
      std::nth_element(sorted.begin(), sorted.begin() + prune_size, sorted.end(), greater_first<value_type>());
      
      const weight_type cutoff = sorted[prune_size].first;
      
      typename sorted_type::const_iterator siter = sorted.begin();
      typename sorted_type::const_iterator siter_end = sorted.end();
      typename sorted_type::const_iterator siter_last = siter + prune_size;

      removed_type removed(source.edges.size(), false);
      
      bool found_equal = false;
      for (/**/; siter != siter_last; ++ siter)
	found_equal |= (siter->first == cutoff);
      
      if (found_equal) {
	for (/**/; siter != siter_end; ++ siter)
	  removed[siter->second] = (siter->first != cutoff);
      } else {
	for (/**/; siter != siter_end; ++ siter)
	  removed[siter->second] = true;
      }
      
      topologically_sort(source, target, filter_pruned(removed), validate);
    }

    const function_type& function;
    const double threshold;
    const bool validate;
  };
  
  
  template <typename Function>
  inline
  void prune_density(const HyperGraph& source, HyperGraph& target, const Function& func, const double threshold, const bool validate=true)
  {
    PruneDensity<Function> __prune(func, threshold, validate);
    
    __prune(source, target);
  }
  
  template <typename Function>
  inline
  void prune_density(HyperGraph& source, const Function& func, const double threshold, const bool validate=true)
  {
    PruneDensity<Function> __prune(func, threshold, validate);

    HyperGraph target;
    
    __prune(source, target);
    
    source.swap(target);
  }
  
};


#endif
