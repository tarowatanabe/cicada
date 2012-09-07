// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__DOT_PRODUCT__HPP__
#define __CICADA__DOT_PRODUCT__HPP__ 1

#include <numeric>

#include <cicada/feature_vector.hpp>
#include <cicada/feature_vector_linear.hpp>
#include <cicada/feature_vector_compact.hpp>
#include <cicada/weight_vector.hpp>

#include <utils/bithack.hpp>

namespace cicada
{
  namespace details
  {
    template <typename Iterator, typename Tp>
    inline
    Tp __inner_product(Iterator first, Iterator last, Tp __dot)
    {
      for (/**/; first != last; ++ first)
	__dot += first->second * first->second;
      return __dot;
    }

    template <typename Iterator, typename Tp1, typename Alloc1, typename Tp>
    inline
    Tp __inner_product(Iterator first, Iterator last, const WeightVector<Tp1, Alloc1>& w, Tp __dot)
    {
      for (/**/; first != last; ++ first)
	__dot += first->second * w[first->first] * first->second;
      return __dot;
    }
  };

  // ordered iterator...
  template <typename Iterator1, typename Iterator2, typename Tp>
  inline
  Tp dot_product(Iterator1 first1, Iterator1 last1, Iterator2 first2, Iterator2 last2, Tp __dot)
  {
    Iterator1 prev1 = last1;
    Iterator2 prev2 = last2;
    while (first1 != last1 && first2 != last2) {
      
      if (prev1 != last1)
	if (prev1->first > first1->first)
	  throw std::runtime_error("not ordered?");

      if (prev2 != last2)
	if (prev2->first > first2->first)
	  throw std::runtime_error("not ordered?");

      if (first1->first < first2->first) {
	prev1 = first1;
        ++ first1;
      } else if (first2->first < first1->first) {
	prev2 = first2;
        ++ first2;
      } else {
        __dot += first1->second * first2->second;
	
	prev1 = first1;
	prev2 = first2;
	
        ++ first1;
        ++ first2;
      }
    }
    return __dot;
  }
  
  template <typename Iterator1, typename Iterator2, typename Tp, typename BinaryOp>
  inline
  Tp dot_product(Iterator1 first1, Iterator1 last1, Iterator2 first2, Iterator2 last2, Tp __dot, BinaryOp op)
  {
    typedef typename std::iterator_traits<Iterator1>::value_type::second_type value1_type;
    typedef typename std::iterator_traits<Iterator2>::value_type::second_type value2_type;
    
    Iterator1 prev1 = last1;
    Iterator2 prev2 = last2;
    
    while (first1 != last1 && first2 != last2) {
      if (prev1 != last1)
	if (prev1->first > first1->first)
	  throw std::runtime_error("not ordered?");

      if (prev2 != last2)
	if (prev2->first > first2->first)
	  throw std::runtime_error("not ordered?");
      
      if (first1->first < first2->first) {
        __dot += op(first1->second, value2_type());
	prev1 = first1;
        ++ first1;
      } else if (first2->first < first1->first) {
        __dot += op(value1_type(), first2->second);
	prev2 = first2;
        ++ first2;
      } else {
        __dot += op(first1->second, first2->second);
        
	prev1 = first1;
	prev2 = first2;

        ++ first1;
        ++ first2;
      } 
    }
    
    for (/**/; first1 != last1; ++ first1)
      __dot += op(first1->second, value2_type());
    for (/**/; first2 != last2; ++ first2)
      __dot += op(value1_type(), first2->second);

    return __dot;
  }
  
  template <typename Iterator1, typename Tp1, typename Alloc1, typename Iterator2, typename Tp>
  inline
  Tp dot_product(Iterator1 first1, Iterator1 last1, const WeightVector<Tp1, Alloc1>& w, Iterator2 first2, Iterator2 last2, Tp __dot)
  {
    Iterator1 prev1 = last1;
    Iterator2 prev2 = last2;

    while (first1 != last1 && first2 != last2) {
      if (prev1 != last1)
	if (prev1->first > first1->first)
	  throw std::runtime_error("not ordered?");

      if (prev2 != last2)
	if (prev2->first > first2->first)
	  throw std::runtime_error("not ordered?");

      if (first1->first < first2->first) {
	prev1 = first1;
        ++ first1;
      } else if (first2->first < first1->first) {
	prev2 = first2;
        ++ first2;
      } else {
        __dot += first1->second * w[first1->first] * first2->second;
        
	prev1 = first1;
	prev2 = first2;

        ++ first1;
        ++ first2;
      } 
    }
    return __dot;
  }
  
  // weight-vector dot weight-vector
  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const WeightVector<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& y)
  {
    const size_t size = utils::bithack::min(x.size(), y.size());
    
    return std::inner_product(x.begin(), x.begin() + size, y.begin(), Tp1());
  }

  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2, typename BinaryOp>
  inline
  Tp1 dot_product(const WeightVector<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& y, BinaryOp op)
  {
    typedef WeightVector<Tp1, Alloc1> weight_vector1_type;
    typedef WeightVector<Tp2, Alloc2> weight_vector2_type;
    
    const size_t size = utils::bithack::min(x.size(), y.size());
    
    Tp1 __dot =  std::inner_product(x.begin(), x.begin() + size, y.begin(), Tp1(), std::plus<Tp1>(), op);
    
    if (x.size() > y.size()) {
      typename weight_vector1_type::const_iterator iter_end = x.end();
      for (typename weight_vector1_type::const_iterator iter = x.begin() + size; iter != iter_end; ++ iter)
	__dot += op(*iter, Tp2());
    } else {
      typename weight_vector2_type::const_iterator iter_end = y.end();
      for (typename weight_vector2_type::const_iterator iter = y.begin() + size; iter != iter_end; ++ iter)
	__dot += op(Tp1(), *iter);
    }
    return __dot;
  }

  template <typename Iterator, typename Tp2, typename Alloc2, typename Tp>
  inline
  Tp dot_product(Iterator first, Iterator last, const WeightVector<Tp2, Alloc2>& y, Tp __dot)
  {
    for (/**/; first != last; ++ first)
      __dot += first->second * y[first->first];
    return __dot;
  }
  
  template <typename Tp1, typename Alloc1, typename Iterator, typename Tp>
  inline
  Tp dot_product(const WeightVector<Tp1, Alloc1>& x, Iterator first, Iterator last, Tp __dot)
  {
    for (/**/; first != last; ++ first)
      __dot += x[first->first] * first->second;
    return __dot;
  }

  // feature-vector dot weight-vector

  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& y)
  {
    return dot_product(x.begin(), x.end(), y, Tp1());
  }
  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const WeightVector<Tp1, Alloc1>& x, const FeatureVector<Tp2, Alloc2>& y)
  {
    return dot_product(x, y.begin(), y.end(), Tp1());
  }
  
  template <typename Iterator, typename Tp2, typename Alloc2, typename Tp>
  inline
  Tp dot_product(Iterator first, Iterator last, const FeatureVector<Tp2, Alloc2>& y, Tp __dot)
  {
    typedef FeatureVector<Tp2, Alloc2> feature_vector2_type;
    
    for (/**/; first != last; ++ first) {
      typename feature_vector2_type::const_iterator iter = y.find(first->first);
      
      if (iter != y.end())
	__dot += first->second * iter->second;
    }
    
    return __dot;
  }

  template <typename Tp1, typename Alloc1, typename Iterator, typename Tp>
  inline
  Tp dot_product(const FeatureVector<Tp1, Alloc1>& x, Iterator first, Iterator last, Tp __dot)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector1_type;
    
    for (/**/; first != last; ++ first) {
      typename feature_vector1_type::const_iterator iter = x.find(first->first);
      
      if (iter != x.end())
	__dot += iter->second * first->second;
    }
    return __dot;
  }

  template <typename Iterator, typename Tp1, typename Alloc1, typename Tp2, typename Alloc2, typename Tp>
  inline
  Tp dot_product(Iterator first, Iterator last, const WeightVector<Tp1, Alloc1>& w, const FeatureVector<Tp2, Alloc2>& y, Tp __dot)
  {
    typedef WeightVector<Tp1, Alloc1> weight_vector_type;
    typedef FeatureVector<Tp2, Alloc2> feature_vector2_type;
    
    for (/**/; first != last; ++ first) {
      typename feature_vector2_type::const_iterator iter = y.find(first->first);
      
      if (iter != y.end())
	__dot += first->second * w[first->first] * iter->second;
    }
    
    return __dot;
  }

  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2, typename Iterator, typename Tp>
  inline
  Tp dot_product(const FeatureVector<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& w, Iterator first, Iterator last, Tp __dot)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector1_type;
    typedef WeightVector<Tp2, Alloc2> weight_vector_type;
    
    for (/**/; first != last; ++ first) {
      typename feature_vector1_type::const_iterator iter = x.find(first->first);
      
      if (iter != x.end())
	__dot += iter->second * w[first->first] * first->second;
    }
    return __dot;
  }
  
  template <typename Tp2, typename Alloc2>
  inline
  FeatureVectorCompact::data_type dot_product(const FeatureVectorCompact& x, const WeightVector<Tp2, Alloc2>& y)
  {
    return dot_product(x.begin(), x.end(), y, FeatureVectorCompact::data_type());
  }
  
  template <typename Tp1, typename Alloc1>
  inline
  Tp1 dot_product(const WeightVector<Tp1, Alloc1>& x, const FeatureVectorCompact& y)
  {
    return dot_product(x, y.begin(), y.end(), Tp1());
  }
  

  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVectorLinear<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& y)
  {
    return dot_product(x.begin(), x.end(), y, Tp1());
  }
  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const WeightVector<Tp1, Alloc1>& x, const FeatureVectorLinear<Tp2, Alloc2>& y)
  {
    return dot_product(x, y.begin(), y.end(), Tp1());
  }

  template <typename Iterator, typename Tp2, typename Alloc2, typename Tp>
  inline
  Tp dot_product(Iterator first, Iterator last, const FeatureVectorLinear<Tp2, Alloc2>& y, Tp __dot)
  {
    typedef FeatureVectorLinear<Tp2, Alloc2> feature_vector2_type;
    
    for (/**/; first != last; ++ first) {
      typename feature_vector2_type::const_iterator iter = y.find(first->first);
      
      if (iter != y.end())
	__dot += first->second * iter->second;
    }
    
    return __dot;
  }

  template <typename Tp1, typename Alloc1, typename Iterator, typename Tp>
  inline
  Tp dot_product(const FeatureVectorLinear<Tp1, Alloc1>& x, Iterator first, Iterator last, Tp __dot)
  {
    typedef FeatureVectorLinear<Tp1, Alloc1> feature_vector1_type;
    
    for (/**/; first != last; ++ first) {
      typename feature_vector1_type::const_iterator iter = x.find(first->first);
      
      if (iter != x.end())
	__dot += iter->second * first->second;
    }
    return __dot;
  }

  template <typename Iterator, typename Tp1, typename Alloc1, typename Tp2, typename Alloc2, typename Tp>
  inline
  Tp dot_product(Iterator first, Iterator last, const WeightVector<Tp1, Alloc1>& w, const FeatureVectorLinear<Tp2, Alloc2>& y, Tp __dot)
  {
    typedef WeightVector<Tp1, Alloc1> weight_vector_type;
    typedef FeatureVectorLinear<Tp2, Alloc2> feature_vector2_type;
    
    for (/**/; first != last; ++ first) {
      typename feature_vector2_type::const_iterator iter = y.find(first->first);
      
      if (iter != y.end())
	__dot += first->second * w[first->first] * iter->second;
    }
    
    return __dot;
  }

  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2, typename Iterator, typename Tp>
  inline
  Tp dot_product(const FeatureVectorLinear<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& w, Iterator first, Iterator last, Tp __dot)
  {
    typedef WeightVector<Tp2, Alloc2> weight_vector_typ;
    typedef FeatureVectorLinear<Tp1, Alloc1> feature_vector1_type;
    
    for (/**/; first != last; ++ first) {
      typename feature_vector1_type::const_iterator iter = x.find(first->first);
      
      if (iter != x.end())
	__dot += iter->second * w[first->first] * first->second;
    }
    return __dot;
  }
  
  // feature-vector dot feature-vector
  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const FeatureVector<Tp2, Alloc2>& y)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector1_type;
    typedef FeatureVector<Tp2, Alloc2> feature_vector2_type;
    
    if (x.empty() || y.empty()) return Tp1();
    
    // if the same address, we are identical!
    if (static_cast<const void*>(&x) == static_cast<const void*>(&y))
      return details::__inner_product(x.begin(), x.end(), Tp1());
    
    if (x.size() < y.size())
      return dot_product(x.begin(), x.end(), y, Tp1());
    else
      return dot_product(x, y.begin(), y.end(), Tp1());
  }
  
  template <typename Tp1, typename Alloc1, typename Tp, typename Alloc, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const WeightVector<Tp, Alloc>& w, const FeatureVector<Tp2, Alloc2>& y)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector1_type;
    typedef WeightVector<Tp, Alloc>    weight_vector_type;
    typedef FeatureVector<Tp2, Alloc2> feature_vector2_type;
    
    if (x.empty() || y.empty()) return Tp1();
    
    if (static_cast<const void*>(&x) == static_cast<const void*>(&y))
      return details::__inner_product(x.begin(), x.end(), w, Tp1());
    
if (x.size() < y.size())
      return dot_product(x.begin(), x.end(), w, y, Tp1());
    else
      return dot_product(x, w, y.begin(), y.end(), Tp1());
  }

};

#endif
