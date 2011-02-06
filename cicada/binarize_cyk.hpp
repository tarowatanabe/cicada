// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_CYK__HPP__
#define __CICADA__BINARIZE_CYK__HPP__ 1

#include <cicada/binarize_base.hpp>

namespace cicada
{
  struct BinarizeCYK : public BinarizeBase
  {
    BinarizeCYK(const int __order=-1)
      : order(__order) {}

    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, copy...
      target = source;
      
      
      
    }
    
    const int order;
  };

  inline
  void binarize_cyk(const HyperGraph& source, HyperGraph& target, const int order=-1)
  {
    BinarizeCYK binarizer(order);
    
    binarizer(source, target);
  }

  inline
  void binarize_cyk(HyperGraph& source, const int order=-1)
  {
    HyperGraph target;

    BinarizeCYK binarizer(order);
    
    binarizer(source, target);
    
    source.swap(target);
  }
};

#endif
