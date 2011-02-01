// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PRUNE_BEAM__HPP__
#define __CICADA__PRUNE_BEAM__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sort.hpp>
#include <cicada/inside_outside.hpp>

namespace cicada
{
  template <typename Function>
  struct PruneBeam
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


    PruneBeam(const function_type& __function,
	      const double __threshold,
	      const bool __validate=true)
      : function(__function),
	threshold(__threshold),
	validate(__validate) {}
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      typedef std::vector<weight_type, std::allocator<weight_type> > inside_type;
      typedef std::vector<weight_type, std::allocator<weight_type> > posterior_type;
      
      target.clear();
      
      if (! source.is_valid())
	return;
      
      inside_type    inside(source.nodes.size());
      posterior_type posterior(source.edges.size());
      
      inside_outside(source, inside, posterior, function, function);
      
      // compute max...
      weight_type posterior_max;
      for (id_type id = 0; id != source.edges.size(); ++ id)
	posterior_max = std::max(posterior_max, posterior[id]);
      
      const weight_type cutoff(posterior_max * cicada::semiring::traits<weight_type>::exp(- threshold));
      
      removed_type removed(source.edges.size(), false);
      for (id_type id = 0; id != source.edges.size(); ++ id)
	removed[id] = (posterior[id] < cutoff);
      
      topologically_sort(source, target, filter_pruned(removed), validate);
    }

    const function_type& function;
    const double threshold;
    const bool validate;
  };
  
  
  template <typename Function>
  inline
  void prune_beam(const HyperGraph& source, HyperGraph& target, const Function& func, const double threshold, const bool validate=true)
  {
    PruneBeam<Function> __prune(func, threshold, validate);
    
    __prune(source, target);
  }
  
  template <typename Function>
  inline
  void prune_beam(HyperGraph& source, const Function& func, const double threshold, const bool validate=true)
  {
    PruneBeam<Function> __prune(func, threshold, validate);

    HyperGraph target;
    
    __prune(source, target);
    
    source.swap(target);
  }

};

#endif
