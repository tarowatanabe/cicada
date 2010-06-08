// -*- mode: c++ -*-

#ifndef __CICADA__PRUNE__HPP__
#define __CICADA__PRUNE__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sort.hpp>
#include <cicada/inside_outside.hpp>

namespace cicada
{
  template <typename WeightSet>
  struct BeamPrune
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    typedef WeightSet weight_set_type;
    typedef cicada::semiring::Logprob<double> weight_type;
    
    struct weight_function
    {
      typedef weight_type value_type;
      
      weight_function(const double __scale,
		      const weight_set_type& __weights)
	: scale(__scale), weights(__weights) {}
      
      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	return weight_type::log(edge.features.dot(weights) * scale);
      }
      
      const double& scale;
      const weight_set_type& weights;
    };
    
    
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


    BeamPrune(const weight_set_type& __weights,
	      const double __scale,
	      const double __threshold)
      : weights(__weights),
	scale(__scale),
	threshold(__threshold) {}
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      typedef std::vector<weight_type, std::allocator<weight_type> > inside_type;
      typedef std::vector<weight_type, std::allocator<weight_type> > posterior_type;

      if (source.goal == hypergraph_type::invalid)
	throw std::runtime_error("invalid graph");

      target.clear();
      
      inside_type    inside(source.nodes.size());
      posterior_type posterior(source.edges.size());
      
      inside_outside(source, inside, posterior, weight_function(scale, weights), weight_function(scale, weights));
      
      // compute max...
      weight_type posterior_max;
      for (id_type id = 0; id != source.edges.size(); ++ id)
	posterior_max = std::max(posterior_max, posterior[id]);
      
      const weight_type cutoff(posterior_max * weight_type(threshold));
      
      removed_type removed(source.edges.size(), false);
      size_t num_removed = 0;
      for (id_type id = 0; id != source.edges.size(); ++ id) {
	removed[id] = (posterior[id] < cutoff);
	num_removed += (posterior[id] < cutoff);
      }
      
      topologically_sort(source, target, filter_pruned(removed));
    }

    const weight_set_type& weights;
    const double scale;
    const double threshold;
  };
  
  
  template <typename WeightSet>
  inline
  void beam_prune(const HyperGraph& source, HyperGraph& target, const WeightSet& weights, const double scale, const double threshold)
  {
    BeamPrune<WeightSet> __prune(weights, scale, threshold);
    
    __prune(source, target);
  }
  
  template <typename WeightSet>
  inline
  void beam_prune(HyperGraph& source, const WeightSet& weights, const double scale, const double threshold)
  {
    BeamPrune<WeightSet> __prune(weights, scale, threshold);

    HyperGraph target;
    
    __prune(source, target);
    
    source.swap(target);
  }

};

#endif
