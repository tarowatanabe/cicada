// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__DOT_PRODUCT__HPP__
#define __CICADA__DOT_PRODUCT__HPP__ 1

#include <cicada/feature_vector.hpp>
#include <cicada/weight_vector.hpp>

namespace cicada
{
  
  template <typename Tp, typename Alloc>
  inline
  Tp dot_product(const FeatureVector<Tp, Alloc>& x)
  {
    typedef FeatureVector<Tp, Alloc> feature_vector_type;
    
    Tp sum = Tp();
    typename feature_vector_type::const_iterator iter_end = x.end();
    for (typename feature_vector_type::const_iterator iter = x.begin(); iter != iter_end; ++ iter)
      sum += iter->second;
    
    return sum;
  }
  
};

#endif
