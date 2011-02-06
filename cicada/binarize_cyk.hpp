// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_CYK__HPP__
#define __CICADA__BINARIZE_CYK__HPP__ 1

#include <deque>

#include <cicada/binarize_base.hpp>

namespace cicada
{
  struct BinarizeCYK : public BinarizeBase
  {
    BinarizeCYK(const int __order=0)
      : order(__order) {}
    
    typedef std::pair<int, int> span_type;
    typedef std::vector<span_type, std::allocator<span_type> > span_set_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > label_set_type;
    typedef std::deque<label_set_type, std::allocator<label_set_type> > label_nodes_type;
    
    typedef utils::chart<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_chart_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, topological-order traversal to compute span for each node and 
      // ancestors.
      // I don't know how to handle unary rules... we will take label from the bottom. (top or concatenate them together?)
      
      // order <= 0 implies all the ancestors...
      
    }
    
    span_set_type spans;
    
    const int order;
  };
  
  inline
  void binarize_cyk(const HyperGraph& source, HyperGraph& target, const int order=0)
  {
    BinarizeCYK binarizer(order);
    
    binarizer(source, target);
  }

  inline
  void binarize_cyk(HyperGraph& source, const int order=0)
  {
    HyperGraph target;

    BinarizeCYK binarizer(order);
    
    binarizer(source, target);
    
    source.swap(target);
  }
};

#endif
