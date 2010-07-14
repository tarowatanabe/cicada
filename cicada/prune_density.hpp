// -*- mode: c++ -*-

#ifndef __CICADA__PRUNE_DENSITY__HPP__
#define __CICADA__PRUNE_DENSITY__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sort.hpp>
#include <cicada/viterbi.hpp>
#include <cicada/inside_outside.hpp>

#include <utils/sgi_hash_set.hpp>

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

    struct edge_traversal
    {
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > value_type;
      
      template <typename Edge, typename Iterator>
      void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
      {
	yield.clear();
	yield.push_back(edge.id);
	for (/**/; first != last; ++ first)
	  yield.insert(yield.end(), first->begin(), first->end());
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

#ifdef HAVE_TR1_UNORDERED_SET
      typedef std::tr1::unordered_set<id_type, utils::hashmurmur<size_t>, std::equal_to<id_type>,
	std::allocator<id_type> > id_set_type;
#else
      typedef sgi::hash_set<id_type, utils::hashmurmur<size_t>, std::equal_to<id_type>,
	std::allocator<id_type> > id_set_type;
#endif

      
      if (source.goal == hypergraph_type::invalid)
	throw std::runtime_error("invalid graph");
      
      target.clear();
      
      weight_type viterbi_weight;
      edge_set_type viterbi_edges;
      viterbi(source, viterbi_edges, viterbi_weight, edge_traversal(), function);
      
      const size_t prune_size = size_t(threshold * viterbi_edges.size());
      
      if (source.edges.size() <= prune_size) {
	target = source;
	return;
      }

      id_set_type edges_unique(viterbi_edges.begin(), viterbi_edges.end());
      
      inside_type    inside(source.nodes.size());
      posterior_type posterior(source.edges.size());
      
      inside_outside(source, inside, posterior, function, function);
      
      posterior_type sorted(posterior);
      
      std::nth_element(sorted.begin(), sorted.begin() + prune_size, sorted.end(), std::greater<weight_type>());
      
      const weight_type cutoff = *(sorted.begin() + prune_size);
      
      removed_type removed(source.edges.size(), false);
      size_t num_removed = 0;
      for (id_type id = 0; id != source.edges.size(); ++ id) {
	const bool remove = (posterior[id] < cutoff) && (edges_unique.find(id) == edges_unique.end());
	
	removed[id] = remove;
	num_removed += remove;
      }
      
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
