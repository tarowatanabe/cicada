// -*- mode: c++ -*-

#ifndef __CICADA__PRUNE_DENSITY__HPP__
#define __CICADA__PRUNE_DENSITY__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sort.hpp>
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
		 const double __threshold)
      : function(__function),
	threshold(__threshold) {}
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      typedef std::vector<weight_type, std::allocator<weight_type> > inside_type;
      typedef std::vector<weight_type, std::allocator<weight_type> > posterior_type;

      typedef hypergraph_type::id_type id_type;
      
      typedef std::vector<id_type, std::allocator<id_type> > edge_set_type;
      
      if (source.goal == hypergraph_type::invalid)
	throw std::runtime_error("invalid graph");
      
      target.clear();
      
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
      
      posterior_type sorted(posterior);
      
      std::nth_element(sorted.begin(), sorted.begin() + prune_size, sorted.end(), std::greater<weight_type>());
      
      // 1e-7 for adjusting numerical instability
      const weight_type cutoff = sorted[prune_size] * cicada::semiring::traits<weight_type>::log(- 1e-7);
      
      removed_type removed(source.edges.size(), false);
      for (id_type id = 0; id != source.edges.size(); ++ id)
	removed[id] = (posterior[id] < cutoff);
      
      topologically_sort(source, target, filter_pruned(removed));
    }

    const function_type& function;
    const double threshold;
  };
  
  
  template <typename Function>
  inline
  void prune_density(const HyperGraph& source, HyperGraph& target, const Function& func, const double threshold)
  {
    PruneDensity<Function> __prune(func, threshold);
    
    __prune(source, target);
  }
  
  template <typename Function>
  inline
  void prune_density(HyperGraph& source, const Function& func, const double threshold)
  {
    PruneDensity<Function> __prune(func, threshold);

    HyperGraph target;
    
    __prune(source, target);
    
    source.swap(target);
  }
  
};


#endif
