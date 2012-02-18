// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PRUNE_KBEST__HPP__
#define __CICADA__PRUNE_KBEST__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sort_topologically.hpp>
#include <cicada/kbest.hpp>
#include <cicada/inside_outside.hpp>

#include <google/dense_hash_set>

namespace cicada
{
  template <typename Function>
  struct PruneKBest
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

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

    struct traversal
    {
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;
      
      typedef edge_set_type value_type;

      template <typename Edge, typename Iterator>
      void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
      {
	yield.clear();
	
	yield.push_back(edge.id);
	for (/**/; first != last; ++ first)
	  yield.insert(yield.end(), first->begin(), first->end());
      }
    };
    
    struct kbest_filter
    {
      template <typename Node, typename Yield>
      bool operator()(const Node& node, const Yield& yield) const
      {
	return false;
      }
    };

    typedef google::dense_hash_set<hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<hypergraph_type::id_type> > edge_set_type;
    
    PruneKBest(const function_type& __function,
	       const size_type __kbest_size,
	       const bool __validate=true)
      : function(__function),
	kbest_size(__kbest_size),
	validate(__validate) {}
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      typedef std::vector<weight_type, std::allocator<weight_type> > inside_type;
      typedef std::vector<weight_type, std::allocator<weight_type> > posterior_type;

      typedef cicada::KBest<traversal, Function, kbest_filter> kbest_derivations_type;

      target.clear();
      if (! source.is_valid())
	return;
      
      kbest_derivations_type derivations(source, kbest_size, traversal(), function, kbest_filter());
      
      edge_set_type edges;
      edges.set_empty_key(hypergraph_type::invalid);
      
      typename traversal::value_type derivation;
      size_type   k = 0;
      for (/**/; k != kbest_size; ++ k) {
	weight_type weight;
	
	if (! derivations(k, derivation, weight))
	  break;
	
	edges.insert(derivation.begin(), derivation.end());
      }
      
      if (k != kbest_size) {
	target = source;
	return;
      }
      
      inside_type    inside(source.nodes.size());
      posterior_type posterior(source.edges.size());
      
      inside_outside(source, inside, posterior, function, function);
      
      weight_type weight(cicada::semiring::traits<weight_type>::max());
      
      typename edge_set_type::const_iterator eiter_end = edges.end();
      for (typename edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter)
	weight = std::min(weight, posterior[*eiter]);
      
      removed_type removed(source.edges.size(), false);
      for (id_type id = 0; id != source.edges.size(); ++ id)
	removed[id] = (posterior[id] < weight);
      
      topologically_sort(source, target, filter_pruned(removed), validate);

      if (! target.is_valid())
	target = source;
    }
    
    const function_type& function;
    const size_type kbest_size;
    const bool validate;
  };
  
  
  template <typename Function>
  inline
  void prune_kbest(const HyperGraph& source, HyperGraph& target, const Function& func, const size_t kbest_size, const bool validate=true)
  {
    PruneKBest<Function> __prune(func, kbest_size, validate);
    
    __prune(source, target);
  }
  
  template <typename Function>
  inline
  void prune_kbest(HyperGraph& source, const Function& func, const size_t kbest_size, const bool validate=true)
  {
    PruneKBest<Function> __prune(func, kbest_size, validate);

    HyperGraph target;
    
    __prune(source, target);
    
    source.swap(target);
  }

};

#endif
