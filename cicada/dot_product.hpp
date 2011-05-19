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

  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& y)
  {
    return (x.sparse()
	    ? dot_product(x.sparse_begin(), x.sparse_end(), y, Tp1())
	    : dot_product(x.dense_begin(), x.dense_end(), y, Tp1()));
  }
  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const WeightVector<Tp1, Alloc1>& x, const FeatureVector<Tp2, Alloc2>& y)
  {
    return (y.sparse()
	    ? dot_product(x, y.sparse_begin(), y.sparse_end(), Tp1())
	    : dot_product(x, y.dense_begin(), y.dense_end(), Tp1()));
  }

  template <typename Iterator1, typename Iterator2, typename Tp>
  inline
  Tp dot_product(Iterator1 first1, Iterator1 last1, Iterator2 first2, Iterator2 last2, Tp __dot)
  {
    while (first1 != last1 && first2 != last2) {
      if (first1->first < first2->first)
	++ first1;
      else if (first2->first < first1->first)
	++ first2;
      else {
	__dot += first1->second * first2->second;
	
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
    
    while (first1 != last1 && first2 != last2) {
      if (first1->first < first2->first) {
	__dot += op(first1->second, value2_type());
	++ first1;
      } else if (first2->first < first1->first) {
	__dot += op(value1_type(), first2->second);
	++ first2;
      } else {
	__dot += op(first1->second, first2->second);
	
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

  template <typename Iterator, typename Tp>
  inline
  Tp __dot_product(Iterator first, Iterator last, Tp __dot)
  {
    for (/**/; first != last; ++ first)
      __dot += first->second * first->second;
    return __dot;
  }

  template <typename Iterator, typename Tp, typename BinaryOp>
  inline
  Tp __dot_product(Iterator first, Iterator last, Tp __dot, BinaryOp op)
  {
    for (/**/; first != last; ++ first)
      __dot += op(first->second, first->second);
    return __dot;
  }
  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const FeatureVector<Tp2, Alloc2>& y)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector1_type;
    typedef FeatureVector<Tp2, Alloc2> feature_vector2_type;
    
    if (x.empty() || y.empty()) return Tp1();
    
    // if the same address, we are identical!
    if (static_cast<const void*>(&x) == static_cast<const void*>(&y))
      return (x.sparse()
	      ? __dot_product(x.sparse_begin(), x.sparse_end(), Tp1())
	      : __dot_product(x.dense_begin(), x.dense_end(), Tp1()));
    
    if (x.sparse()) {
      typename feature_vector1_type::const_sparse_iterator iter1     = x.lower_bound(y.begin()->first);
      typename feature_vector1_type::const_sparse_iterator iter1_end = x.end();
      
      if (y.sparse()) {
	typename feature_vector2_type::const_sparse_iterator iter2 = (iter1 != iter1_end ? y.lower_bound(iter1->first) : y.end());
	typename feature_vector2_type::const_sparse_iterator iter2_end = y.end();
	
	return dot_product(iter1, iter1_end, iter2, iter2_end, Tp1());
      } else {
	typename feature_vector2_type::const_dense_iterator iter2 = (iter1 != iter1_end ? y.lower_bound(iter1->first) : y.end());
	typename feature_vector2_type::const_dense_iterator iter2_end = y.end();
	
	return dot_product(iter1, iter1_end, iter2, iter2_end, Tp1());
      }
    } else {
      typename feature_vector1_type::const_dense_iterator iter1     = x.lower_bound(y.begin()->first);
      typename feature_vector1_type::const_dense_iterator iter1_end = x.end();
      
      if (y.sparse()) {
	typename feature_vector2_type::const_sparse_iterator iter2 = (iter1 != iter1_end ? y.lower_bound(iter1->first) : y.end());
	typename feature_vector2_type::const_sparse_iterator iter2_end = y.end();
	
	return dot_product(iter1, iter1_end, iter2, iter2_end, Tp1());
      } else {
	typename feature_vector2_type::const_dense_iterator iter2 = (iter1 != iter1_end ? y.lower_bound(iter1->first) : y.end());
	typename feature_vector2_type::const_dense_iterator iter2_end = y.end();
	
	return dot_product(iter1, iter1_end, iter2, iter2_end, Tp1());
      }
    }
  }

  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2, typename BinaryOp>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const FeatureVector<Tp2, Alloc2>& y, BinaryOp op)
  {
    if (x.empty() && y.empty()) return Tp1();
    
    if (static_cast<const void*>(&x) == static_cast<const void*>(&y))
      return (x.sparse()
	      ? __dot_product(x.sparse_begin(), x.sparse_end(), Tp1(), op)
	      : __dot_product(x.dense_begin(), x.dense_end(), Tp1(), op));
    
    if (x.sparse())
      return (y.sparse()
	      ? dot_product(x.sparse_begin(), x.sparse_end(), y.sparse_begin(), y.sparse_end(), Tp1(), op)
	      : dot_product(x.sparse_begin(), x.sparse_end(), y.dense_begin(), y.dense_end(), Tp1(), op));
    else
      return (y.sparse()
	      ? dot_product(x.dense_begin(), x.dense_end(), y.sparse_begin(), y.sparse_end(), Tp1(), op)
	      : dot_product(x.dense_begin(), x.dense_end(), y.dense_begin(), y.dense_end(), Tp1(), op));
  }
  
  template <typename Iterator, typename Tp1, typename Alloc1, typename Tp>
  inline
  Tp __dot_product(Iterator first, Iterator last, const WeightVector<Tp1, Alloc1>& w, Tp __dot)
  {
    for (/**/; first != last; ++ first)
      __dot += first->second * w[first->first] * first->second;
    return __dot;
  }
  
  template <typename Iterator1, typename Tp1, typename Alloc1, typename Iterator2, typename Tp>
  inline
  Tp dot_product(Iterator1 first1, Iterator1 last1, const WeightVector<Tp1, Alloc1>& w, Iterator2 first2, Iterator2 last2, Tp __dot)
  {
        while (first1 != last1 && first2 != last2) {
      if (first1->first < first2->first)
	++ first1;
      else if (first2->first < first1->first)
	++ first2;
      else {
	__dot += first1->second * w[first1->first] * first2->second;
	
	++ first1;
	++ first2;
      }
    }
    return __dot;
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
      return (x.sparse()
	      ? __dot_product(x.sparse_begin(), x.sparse_end(), w, Tp1())
	      : __dot_product(x.dense_begin(),  x.dense_end(), w, Tp1()));
    
    if (x.sparse()) {
      typename feature_vector1_type::const_sparse_iterator iter1     = x.lower_bound(y.begin()->first);
      typename feature_vector1_type::const_sparse_iterator iter1_end = x.end();
      
      if (y.sparse()) {
	typename feature_vector2_type::const_sparse_iterator iter2 = (iter1 != iter1_end ? y.lower_bound(iter1->first) : y.end());
	typename feature_vector2_type::const_sparse_iterator iter2_end = y.end();
	
	return dot_product(iter1, iter1_end, w, iter2, iter2_end, Tp1());
      } else {
	typename feature_vector2_type::const_dense_iterator iter2 = (iter1 != iter1_end ? y.lower_bound(iter1->first) : y.end());
	typename feature_vector2_type::const_dense_iterator iter2_end = y.end();
	
	return dot_product(iter1, iter1_end, w, iter2, iter2_end, Tp1());
      }
    } else {
      typename feature_vector1_type::const_dense_iterator iter1     = x.lower_bound(y.begin()->first);
      typename feature_vector1_type::const_dense_iterator iter1_end = x.end();
      
      if (y.sparse()) {
	typename feature_vector2_type::const_sparse_iterator iter2 = (iter1 != iter1_end ? y.lower_bound(iter1->first) : y.end());
	typename feature_vector2_type::const_sparse_iterator iter2_end = y.end();
	
	return dot_product(iter1, iter1_end, w, iter2, iter2_end, Tp1());
      } else {
	typename feature_vector2_type::const_dense_iterator iter2 = (iter1 != iter1_end ? y.lower_bound(iter1->first) : y.end());
	typename feature_vector2_type::const_dense_iterator iter2_end = y.end();
	
	return dot_product(iter1, iter1_end, w, iter2, iter2_end, Tp1());
      }
    }
  }
  

};

#endif
