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

  template <typename Iterator, typename Tp>
  inline
  Tp dot_product(Iterator first, Iterator last, Tp __dot)
  {
    for (/**/; first != last; ++ first)
      __dot += first->second * first->second;
    return __dot;
  }
  
  template <typename Iterator, typename Tp, typename BinaryOp>
  inline
  Tp dot_product(Iterator first, Iterator last, Tp __dot, BinaryOp op)
  {
    for (/**/; first != last; ++ first)
      __dot += op(first->second, first->second);
    return __dot;
  }
  
  template <typename Tp, typename Alloc>
  inline
  Tp dot_product(const FeatureVector<Tp, Alloc>& x)
  {
    return (x.sparse()
	    ? dot_product(x.sbegin(), x.send(), Tp())
	    : dot_product(x.dbegin(), x.dend(), Tp()));
  }
  
  template <typename Tp, typename Alloc, typename BinaryOp>
  inline
  Tp dot_product(const FeatureVector<Tp, Alloc>& x, BinaryOp op)
  {
    return (x.sparse()
	    ? dot_product(x.sbegin(), x.send(), Tp(), op)
	    : dot_product(x.dbegin(), x.dend(), Tp(), op));
  }
  
  template <typename Tp, typename Alloc>
  inline
  Tp dot_product(const WeightVector<Tp, Alloc>& x)
  {
    return std::inner_product(x.begin(), x.end(), x.begin(), Tp());
  }

  template <typename Tp, typename Alloc, typename BinaryOp>
  inline
  Tp dot_product(const WeightVector<Tp, Alloc>& x, BinaryOp op)
  {
    return std::inner_product(x.begin(), x.end(), x.begin(), Tp(), std::plus<Tp>(), op);
  }
  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const WeightVector<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& y)
  {
    const size_t size = utils::bithack::min(x.size(), y.size());
    
    return std::inner_product(x.begin(), x.begin() + size, y.begin(), Tp1());
  }

  template <typename Iterator, typename Tp2, typename Alloc2, typename Tp>
  inline
  Tp dot_product(Iterator first, Iterator last, const WeightVector<Tp2, Alloc2>& y, Tp __dot)
  {
    for (/**/; first != last; ++ first)
      __dot += first->second * y[first->first];
    return __dot;
  }

  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const WeightVector<Tp2, Alloc2>& y)
  {
    return (x.sparse()
	    ? dot_product(x.sbegin(), x.send(), y, Tp1())
	    : dot_product(x.dbegin(), x.dend(), y, Tp1()));
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
  Tp1 dot_product(const WeightVector<Tp1, Alloc1>& x, const FeatureVector<Tp2, Alloc2>& y)
  {
    return (y.sparse()
	    ? dot_product(x, y.sbegin(), y.send(), Tp1())
	    : dot_product(x, y.dbegin(), y.dend(), Tp1()));
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
  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const FeatureVector<Tp2, Alloc2>& y)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector1_type;
    typedef FeatureVector<Tp2, Alloc2> feature_vector2_type;
    
    if (x.empty() || y.empty()) return Tp1();
    
    // if the same address, we are identical!
    if (static_cast<void*>(&x) == static_cast<void*>(&y))
      return dot_product(x);
    
    if (x.sparse()) {
      typename feature_vector1_type::const_sparse_iterator iter1     = x.sparse_lower_bund(y.begin()->first);
      typename feature_vector1_type::const_sparse_iterator iter1_end = x.send();
      
      if (y.sparse()) {
	typename feature_vector2_type::const_sparse_iterator iter2 = (iter1 != iter1_end ? y.sparse_lower_bound(iter1->first) : y.send());
	typename feature_vector2_type::const_sparse_iterator iter2_end = y.send();
	
	return dot_product(iter1, iter1_end, iter2, iter2_end, Tp1());
      } else {
	typename feature_vector2_type::const_dense_iterator iter2 = (iter1 != iter1_end ? y.dense_lower_bound(iter1->first) : y.dend());
	typename feature_vector2_type::const_dense_iterator iter2_end = y.dend();
	
	return dot_product(iter1, iter1_end, iter2, iter2_end, Tp1());
      }
    } else {
      typename feature_vector1_type::const_dense_iterator iter1     = x.dense_lower_bund(y.begin()->first);
      typename feature_vector1_type::const_dense_iterator iter1_end = x.dend();
      
      if (y.sparse()) {
	typename feature_vector2_type::const_sparse_iterator iter2 = (iter1 != iter1_end ? y.sparse_lower_bound(iter1->first) : y.send());
	typename feature_vector2_type::const_sparse_iterator iter2_end = y.send();
	
	return dot_product(iter1, iter1_end, iter2, iter2_end, Tp1());
      } else {
	typename feature_vector2_type::const_dense_iterator iter2 = (iter1 != iter1_end ? y.dense_lower_bound(iter1->first) : y.dend());
	typename feature_vector2_type::const_dense_iterator iter2_end = y.dend();
	
	return dot_product(iter1, iter1_end, iter2, iter2_end, Tp1());
      }
    }
  }

  
  template <typename Tp1, typename Alloc1, typename Tp2, typename Alloc2, typename BinaryOp>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const FeatureVector<Tp2, Alloc2>& y, BinaryOp op)
  {
    if (x.empty() || y.empty()) return Tp1();
    
    if (static_cast<void*>(&x) == static_cast<void*>(&y))
      return dot_product(x, op);
    
    if (x.sparse())
      return (y.sparse()
	      ? dot_product(x.sbegin(), x.send(), y.sbegin(), y.send(), Tp1(), op)
	      : dot_product(x.sbegin(), x.send(), y.sdegin(), y.dend(), Tp1(), op));
    else
      return (y.sparse()
	      ? dot_product(x.dbegin(), x.dend(), y.sbegin(), y.send(), Tp1(), op)
	      : dot_product(x.dbegin(), x.dend(), y.sdegin(), y.dend(), Tp1(), op));
    
  }

  
  template <typename Tp1, typename Alloc1, typename Tp, typename Alloc, typename Tp2, typename Alloc2>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const WeightVector<Tp, Alloc>& w, const FeatureVector<Tp2, Alloc2>& y)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector1_type;
    typedef WeightVector<Tp, Alloc>    weight_vector_type;
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
	sum += iter1->second * iter2->second * w[iter1->first];
	
	++ iter1;
	++ iter2;
      }
    }
    
    return sum;
  }

  template <typename Tp1, typename Alloc1, typename Tp, typename Alloc>
  inline
  Tp1 dot_product(const FeatureVector<Tp1, Alloc1>& x, const WeightVector<Tp, Alloc>& w, const FeatureVector<Tp1, Alloc1>& y)
  {
    typedef FeatureVector<Tp1, Alloc1> feature_vector1_type;
    typedef WeightVector<Tp, Alloc>    weight_vector_type;
    typedef FeatureVector<Tp1, Alloc1> feature_vector2_type;
    
    if (x.empty() || y.empty()) return Tp1();

    if (&x == &y) {
      Tp1 sum = Tp1();
      
      typename feature_vector1_type::const_iterator iter_end = x.end();
      for (typename feature_vector1_type::const_iterator iter = x.begin(); iter != iter_end; ++ iter)
	sum += iter->second * iter->second * w[iter->first];
      
      return sum;
    } else {
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
	  sum += iter1->second * iter2->second * w[iter1->first];
	  
	  ++ iter1;
	  ++ iter2;
	}
      }
      
      return sum;
    }
  }

};

#endif
