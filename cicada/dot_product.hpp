// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__DOT_PRODUCT__HPP__
#define __CICADA__DOT_PRODUCT__HPP__ 1

#include <numeric>

#include <cicada/feature_vector.hpp>
#include <cicada/weight_vector.hpp>

#include <utils/bithack.hpp>

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
      sum += iter->second * iter->second;
    
    return sum;
  }
  
  template <typename Tp, typename Alloc>
  inline
  Tp dot_product(const WeightVector<Tp, Alloc>& x)
  {
    return std::inner_product(x.begin(), x.end(), x.begin(), Tp());
  }
  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const WeightVector<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& y)
  {
    const size_t size = utils::bithack::min(x.size(), y.size());
    
    return std::inner_product(x.begin(), x.begin() + size, y.begin(), Tp1());
  }

  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& y)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector_type;
    
    Tp1 sum = Tp1();
    typename feature_vector_type::const_iterator iter_end = x.end();
    for (typename feature_vector_type::const_iterator iter = x.begin(); iter != iter_end; ++ iter)
      sum += iter->second * y[iter->first];
    
    return sum;
  }
  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const FeatureVector<Tp2, Alloc2>& y)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector1_type;
    typedef FeatureVector<Tp2, Alloc2> feature_vector2_type;
    
    if (x.empty() || y.empty()) return Tp1();

    typename feature_vector1_type::const_iterator iter1     = x.lower_bound(y.begin()->first);
    typename feature_vector1_type::const_iterator iter1_end = x.end();
    
    typename feature_vector2_type::const_iterator iter2     = (iter1 != iter1_end ? y.lower_bound(iter1->first) : y.end());
    typename feature_vector2_type::const_iterator iter2_end = y.end();
    
    Tp1 sum = Tp1();
      
    while (iter1 != iter1_end && iter2 != iter2_end) {
      if (iter1->first < iter2->first)
	++ iter1;
      else if (iter2->first < iter1->first)
	++ iter2;
      else {
	sum += iter1->second * iter2->second;
	
	++ iter1;
	++ iter2;
      }
    }
    
    return sum;
  }

  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2, typename BinaryOp>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const FeatureVector<Tp2, Alloc2>& y, BinaryOp op)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector1_type;
    typedef FeatureVector<Tp2, Alloc2> feature_vector2_type;
    
    if (x.empty() || y.empty()) return Tp1();

    typename feature_vector1_type::const_iterator iter1     = x.begin();
    typename feature_vector1_type::const_iterator iter1_end = x.end();
    
    typename feature_vector2_type::const_iterator iter2     = y.begin();
    typename feature_vector2_type::const_iterator iter2_end = y.end();
    
    Tp1 sum = Tp1();
      
    while (iter1 != iter1_end && iter2 != iter2_end) {
      if (iter1->first < iter2->first) {
	sum += op(iter1->second, Tp2());
	++ iter1;
      } else if (iter2->first < iter1->first) {
	sum += op(Tp1(), iter2->second);
	++ iter2;
      } else {
	sum += op(iter1->second, iter2->second);
	
	++ iter1;
	++ iter2;
      }
    }

    for (/**/; iter1 != iter1_end; ++ iter1)
      sum += op(iter1->second, Tp2());
    for (/**/; iter2 != iter2_end; ++ iter2)
      sum += op(Tp1(), iter2->second);
    
    return sum;
  }

};

#endif
