// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PRUNE_KBEST__HPP__
#define __CICADA__PRUNE_KBEST__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sort.hpp>
#include <cicada/kbest.hpp>

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
    
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;
    typedef std::vector<bool, std::allocator<bool> > survived_type;
    
    struct filter_survived
    {
      const survived_type& survived;
      
      filter_survived(const survived_type& __survived) : survived(__survived) {}
      
      template <typename Edge>
      bool operator()(const Edge& edge) const
      {
	return ! survived[edge.id];
      }
    };

    struct traversal
    {
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
    
    PruneKBest(const function_type& __function,
	       const size_type __kbest_size,
	       const bool __validate=true)
      : function(__function),
	kbest_size(__kbest_size),
	validate(__validate) {}
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      typedef cicada::KBest<traversal, Function, kbest_filter> kbest_derivations_type;

      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      kbest_derivations_type derivations(source, kbest_size, traversal(), function, kbest_filter());
      
      edge_set_type derivation;
      weight_type   weight;
      
      survived_type survived(source.edges.size(), false);
      
      for (int k = 0; k < kbest_size; ++ k) {
	if (! derivations(k, derivation, weight))
	  break;
	
	typename edge_set_type::const_iterator eiter_end = derivation.end();
	for (typename edge_set_type::const_iterator eiter = derivation.begin(); eiter != eiter_end; ++ eiter)
	  survived[*eiter] = true;
      }
      
      cicada::topologically_sort(source, target, filter_survived(survived), validate);
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
