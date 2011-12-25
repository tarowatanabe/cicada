// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE_VECTOR__HPP__
#define __CICADA__FEATURE_VECTOR__HPP__ 1

#include <map>
#include <memory>
#include <utility>
#include <algorithm>
#include <iterator>
#include <iostream>

#include <cicada/feature.hpp>

#include <utils/vector_map.hpp>

namespace cicada
{
  
  // forward declaration...
  template <typename Tp, typename Alloc >
  class WeightVector;

  class FeatureVectorCompact;

  template <typename Tp, typename Alloc >
  class FeatureVectorLinear;

  template <typename Tp, typename Alloc >
  class FeatureVectorUnordered;

  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class FeatureVector
  {
  public:
    typedef cicada::Feature feature_type;
    typedef cicada::Feature key_type;
    typedef Tp mapped_type;
    typedef Tp data_type;
    
    typedef std::pair<const feature_type, data_type> value_type;
    
  private:
    typedef typename Alloc::template rebind<value_type>::other sparse_alloc_type;
    typedef std::map<key_type, data_type, std::less<key_type>, sparse_alloc_type> sparse_vector_type;
    
    typedef std::pair<feature_type, data_type> dense_value_type;
    typedef typename Alloc::template rebind<dense_value_type>::other dense_alloc_type;
    typedef utils::vector_map<key_type, data_type, std::less<key_type>,  dense_alloc_type> dense_vector_type;
    
    typedef FeatureVector<Tp, Alloc> self_type;
    
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef typename sparse_vector_type::const_iterator const_sparse_iterator;
    typedef typename sparse_vector_type::iterator             sparse_iterator;

    typedef typename dense_vector_type::const_iterator const_dense_iterator;
    typedef typename dense_vector_type::iterator             dense_iterator;
    
    template <typename DIterator, typename SIterator, typename Ref, typename Ptr>
    struct __iterator
    {
      typedef std::bidirectional_iterator_tag iterator_category;
      
      typedef Ref       reference;
      typedef Ptr       pointer;
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef std::pair<const feature_type, data_type> value_type;

      __iterator() : diter(0), siter() {}
      __iterator(const typename dense_vector_type::const_iterator& __diter) : diter(__diter), siter() {}
      __iterator(const typename dense_vector_type::iterator& __diter) : diter(__diter), siter() {}
      __iterator(const typename sparse_vector_type::const_iterator& __siter) : diter(0), siter(__siter) {}
      __iterator(const typename sparse_vector_type::iterator& __siter) : diter(0), siter(__siter) {}
      
      template <typename D, typename S, typename R, typename P>
      __iterator(const __iterator<D,S,R,P>& x) : diter(x.diter), siter(x.siter) {}
      __iterator(const __iterator<DIterator,SIterator,Ref,Ptr>& x) : diter(x.diter), siter(x.siter) {}
      
      operator const DIterator() const { return diter; }
      operator const SIterator() const { return siter; }

      template <typename D, typename S, typename R, typename P>
      __iterator& operator=(const __iterator<D,S,R,P>& x)
      {
	diter = x.diter;
	siter = x.siter;
	return *this;
      }
      
      reference operator*() const
      {
	return *(diter ? ((pointer) &(*diter)) : &(*siter));
      }
      
      pointer operator->() const
      {
	return (diter ? ((pointer) &(*diter)) : &(*siter));
      }

      __iterator& operator++()
      {
	if (diter)
	  ++ diter;
	else
	  ++ siter;
	return *this;
      }
      
      __iterator& operator--()
      {
	if (diter)
	  -- diter;
	else
	  -- siter;
	return *this;
      }
      
      __iterator operator++(int)
      {
	__iterator tmp = *this;
	++ *this;
	return tmp;
      }

      __iterator operator--(int)
      {
	__iterator tmp = *this;
	-- *this;
	return tmp;
      }
      
      template <typename D, typename S, typename R, typename P>
      friend
      bool operator==(const __iterator<DIterator,SIterator,Ref,Ptr>& x, const __iterator<D,S,R,P>& y)
      {
	return x.diter == y.diter && x.siter == y.siter;
      }
      
      template <typename D, typename S, typename R, typename P>
      friend
      bool operator!=(const __iterator<DIterator,SIterator,Ref,Ptr>& x, const __iterator<D,S,R,P>& y)
      {
	return x.diter != y.diter || x.siter != y.siter;
      }
      
      DIterator diter;
      SIterator siter;
    };
    
    typedef __iterator<typename dense_vector_type::const_iterator,
		       typename sparse_vector_type::const_iterator,
		       const value_type&,
		       const value_type*> const_iterator;
    typedef __iterator<typename dense_vector_type::iterator,
		       typename sparse_vector_type::iterator,
		       value_type&,
		       value_type*> iterator;
    
    typedef const value_type& const_reference;
    typedef       value_type& reference;
    typedef       value_type* pointer;
    
  private:
    static const size_type __dense_size = 20;
    
  public:
    FeatureVector()
      : __dense(), __sparse(0) {}
    FeatureVector(const FeatureVector<Tp,Alloc>& x)
      : __dense(x.__dense), __sparse(x.__sparse ? new sparse_vector_type(*x.__sparse) : 0) {}
    template <typename T, typename A>
    FeatureVector(const FeatureVector<T,A>& x)
      : __dense(x.__dense.begin(), x.__dense.end()), __sparse(x.__sparse ? new sparse_vector_type(x.__sparse->begin(), x.__sparse->end()) : 0) { }
    template <typename Iterator>
    FeatureVector(Iterator first, Iterator last)
      : __dense(), __sparse(0) { assign(first, last); }
    FeatureVector(const FeatureVectorCompact& x)
      : __dense(), __sparse(0) { assign(x); }
    template <typename T, typename A>
    FeatureVector(const FeatureVectorLinear<T,A>& x)
      : __dense(), __sparse(0) { assign(x); }
    template <typename T, typename A>
    FeatureVector(const FeatureVectorUnordered<T,A>& x)
      : __dense(), __sparse(0) { assign(x); }
    
    FeatureVector& operator=(const FeatureVector<Tp,Alloc>& x)
    {
      assign(x);
      return *this;
    }
    
    template <typename T, typename A>
    FeatureVector& operator=(const FeatureVector<T, A>& x)
    {
      assign(x);
      return *this;
    }

    FeatureVector& operator=(const FeatureVectorCompact& x)
    {
      assign(x);
      return *this;
    }

    template <typename T, typename A>
    FeatureVector& operator=(const FeatureVectorLinear<T, A>& x)
    {
      assign(x);
      return *this;
    }

    template <typename T, typename A>
    FeatureVector& operator=(const FeatureVectorUnordered<T, A>& x)
    {
      assign(x);
      return *this;
    }

    ~FeatureVector()
    {
      if (__sparse) {
	delete __sparse;
	__sparse = 0;
      }
    }
    
  public:
    void assign(const FeatureVector<Tp,Alloc>& x)
    {
      __dense = x.__dense;
      if (x.__sparse) {
	if (__sparse)
	  *__sparse = *x.__sparse;
	else
	  __sparse = new sparse_vector_type(*x.__sparse);
      } else if (__sparse) {
	delete __sparse;
	__sparse = 0;
      }
    }

    template <typename T, typename A>
    void assign(const FeatureVector<T,A>& x)
    {
      __dense.clear();
      __dense.insert(x.__dense.begin(), x.__dense.end());
      if (x.__sparse) {
	if (__sparse) {
	  __sparse->clear();
	  __sparse->insert(x.__sparse->begin(), x.__sparse->end());
	} else
	  __sparse = new sparse_vector_type(x.__sparse->begin(), x.__sparse->end());
      } else if (__sparse) {
	delete __sparse;
	__sparse = 0;
      }
    }

    void assign(const FeatureVectorCompact& x);

    template <typename T, typename A>
    void assign(const FeatureVectorLinear<T,A>& x)
    {
      const size_type __n = x.size();
      
      if (__n > __dense_size) {
	__dense.clear();
	if (__sparse) {
	  __sparse->clear();
	  __sparse->insert(x.begin(), x.end());
	} else
	  __sparse = new sparse_vector_type(x.begin(), x.end());
      } else {
	if (__sparse) {
	  delete __sparse;
	  __sparse = 0;
	}
	__dense.clear();
	__dense.insert(x.begin(), x.end());
      }
    }

    template <typename T, typename A>
    void assign(const FeatureVectorUnordered<T,A>& x)
    {
      const size_type __n = x.size();
      
      if (__n > __dense_size) {
	__dense.clear();
	if (__sparse) {
	  __sparse->clear();
	  __sparse->insert(x.begin(), x.end());
	} else
	  __sparse = new sparse_vector_type(x.begin(), x.end());
      } else {
	if (__sparse) {
	  delete __sparse;
	  __sparse = 0;
	}
	__dense.clear();
	__dense.insert(x.begin(), x.end());
      }
    }
    
    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      const size_type __n = std::distance(first, last);
      
      if (__n > __dense_size) {
	__dense.clear();
	if (__sparse) {
	  __sparse->clear();
	  __sparse->insert(first, last);
	} else
	  __sparse = new sparse_vector_type(first, last);
      } else {
	if (__sparse) {
	  delete __sparse;
	  __sparse = 0;
	}
	__dense.clear();
	__dense.insert(first, last);
      }
    }
    
    template <typename Iterator>
    void insert(Iterator first, Iterator last)
    {
      if (first == last) return;

      if (__sparse)
	__sparse->insert(first, last);
      else {
	__dense.insert(first, last);
	
	if (__dense.size() > __dense_size) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
      }
    }
    
    void insert(iterator iter, const value_type& x)
    {
      if (x.second == Tp()) return;
      
      if (__sparse)
	__sparse->insert(iter.siter, x);
      else {
	__dense.insert(iter.diter, x);
	
	if (__dense.size() > __dense_size) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
      }
    }
    
    void insert(const value_type& x)
    {
      if (x.second == Tp()) return;
      
      if (__sparse)
	__sparse->insert(x);
      else {
	__dense.insert(x);
	
	if (__dense.size() > __dense_size) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
      }
    }

    template <typename T, typename A, typename Prefix>
    void update(const FeatureVector<T,A>& x, const Prefix& prefix)
    {
      if (x.sparse())
	update_ordered(x.sparse_begin(), x.sparse_end(), prefix);
      else
	update_ordered(x.dense_begin(), x.dense_end(), prefix);
    }

    template <typename T, typename A, typename Prefix>
    void update(const FeatureVectorLinear<T,A>& x, const Prefix& prefix)
    {
      update_ordered(x.begin(), x.end(), prefix);
    }

    template <typename T, typename A, typename Prefix>
    void update(const FeatureVectorUnordered<T,A>& x, const Prefix& prefix)
    {
      update_unordered(x.begin(), x.end(), prefix);
    }

  private:
    template <typename Iterator, typename Prefix>
    void update_ordered(Iterator iter2, Iterator iter2_end, const Prefix& prefix)
    {
      if (iter2 == iter2_end)
	erase_prefix(prefix);
      else {
	if (sparse())
	  update_ordered(sparse_begin(), sparse_end(), iter2, iter2_end, prefix);
	else
	  update_ordered(dense_begin(), dense_end(), iter2, iter2_end, prefix);
      }
    }
    
    template <typename Iterator1, typename Iterator2, typename Prefix>
    void update_ordered(Iterator1 iter1, Iterator1 iter1_end, Iterator2 iter2, Iterator2 iter2_end, const Prefix& prefix)
    {
      typedef std::pair<feature_type, data_type> feat_type;
      typedef std::vector<feat_type, std::allocator<feat_type> > feat_set_type;
      
      // update this by x
      // logically, we erase-prefix, then *this += x;
      
      feat_set_type feats;
      feats.reserve(size());
	
      while (iter1 != iter1_end && iter2 != iter2_end) {
	if (iter1->first < iter2->first) {
	  if (iter1->first.size() < prefix.size() || ! std::equal(prefix.begin(), prefix.end(), iter1->first.begin()))
	    feats.push_back(*iter1);
	  ++ iter1;
	} else if (iter2->first < iter1->first) {
	  feats.push_back(*iter2);
	  ++ iter2;
	} else {
	  feats.push_back(*iter2);
	  ++ iter1;
	  ++ iter2;
	}
      }
	
      for (/**/; iter1 != iter1_end; ++ iter1) 
	if (iter1->first.size() < prefix.size() || ! std::equal(prefix.begin(), prefix.end(), iter1->first.begin()))
	  feats.push_back(*iter1);
	
      feats.insert(feats.end(), iter2, iter2_end);
	
      assign(feats.begin(), feats.end());
    }

    template <typename Iterator, typename Prefix>
    void update_unordered(Iterator first, Iterator last, const Prefix& prefix)
    {
      // update this by x
      // logically, we erase-prefix, then *this += x;
      
      erase_prefix(prefix);
      for (/**/; first != last; ++ first)
	operator[](first->first) = first->second;
    }

  public:

    size_type size() const { return (__sparse ? __sparse->size() : __dense.size()); }
    bool empty() const { return (__sparse ? __sparse->empty() : __dense.empty()); }
    bool sparse() const { return __sparse; }

    void reserve(size_type x) { }
    
    void clear()
    {
      __dense.clear();
      
      if (__sparse) {
	delete __sparse;
	__sparse = 0;
      }
    }
    
    Tp operator[](const key_type& x) const
    {
      if (__sparse) {
	typename sparse_vector_type::const_iterator iter = __sparse->find(x);
	return (iter == __sparse->end() ? Tp() : iter->second);
      } else {
	typename dense_vector_type::const_iterator iter = __dense.find(x);
	return (iter == __dense.end() ? Tp() : iter->second);
      }
    }
    
    Tp& operator[](const key_type& x)
    {
      if (__sparse)
	return __sparse->operator[](x);
      else {
	std::pair<typename dense_vector_type::iterator, bool> result = __dense.insert(std::make_pair(x, Tp()));
	if (! result.second || __dense.size() <= __dense_size)
	  return result.first->second;
	
	__sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	__dense.clear();
	return __sparse->operator[](x);
      }
    }
    
    const_reference front() const { return (__sparse ? *__sparse->begin() : *((const value_type*) &__dense.front())); }
    reference front() { return (__sparse ? *__sparse->begin() : *((value_type*) &__dense.front())); }
    
    const_reference back() const { return (__sparse ? *__sparse->rbegin() : *((const value_type*) &__dense.back())); }
    reference back() { return (__sparse ? *__sparse->rbegin() : *((value_type*) &__dense.back()));}

    
    const_iterator find(const key_type& x) const
    {
      return (__sparse ? const_iterator(__sparse->find(x)) : const_iterator(__dense.find(x)));
    }
    
    iterator find(const key_type& x)
    {
      return (__sparse ? iterator(__sparse->find(x)) : iterator(__dense.find(x)));
    }
    
    const_iterator lower_bound(const key_type& x) const
    {
      return (__sparse ? const_iterator(__sparse->lower_bound(x)) : const_iterator(__dense.lower_bound(x)));
    }
    
    iterator lower_bound(const key_type& x)
    {
      return (__sparse ? iterator(__sparse->lower_bound(x)) : iterator(__dense.lower_bound(x)));
    }

    const_iterator upper_bound(const key_type& x) const
    {
      return (__sparse ? const_iterator(__sparse->upper_bound(x)) : const_iterator(__dense.upper_bound(x)));
    }
    
    iterator upper_bound(const key_type& x)
    {
      return (__sparse ? iterator(__sparse->upper_bound(x)) : iterator(__dense.upper_bound(x)));
    }
    
    
    void erase(const key_type& x)
    {
      if (__sparse)
	__sparse->erase(x);
      else
	__dense.erase(x);
    }

    void erase(iterator x)
    {
      if (__sparse)
	__sparse->erase(x.siter);
      else
	__dense.erase(x.diter);
    }
    
    template <typename Prefix>
    void erase_prefix(const Prefix& prefix)
    {
      if (__sparse) {
	for (typename sparse_vector_type::iterator fiter = __sparse->begin(); fiter != __sparse->end(); /**/)
	  if (fiter->first.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), fiter->first.begin()))
	    __sparse->erase(fiter ++);
	  else
	    ++ fiter;
      } else {
	for (typename dense_vector_type::iterator fiter = __dense.begin(); fiter != __dense.end(); ++ fiter)
	  if (fiter->first.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), fiter->first.begin()))
	    fiter->second = Tp();
      }
    }
    
    inline const_iterator begin() const { return (__sparse ? const_iterator(__sparse->begin()) : const_iterator(__dense.begin())); }
    inline       iterator begin()       { return (__sparse ? iterator(__sparse->begin()) : iterator(__dense.begin())); }
    
    inline const_iterator end() const { return (__sparse ? const_iterator(__sparse->end()) : const_iterator(__dense.end())); }
    inline       iterator end()       { return (__sparse ? iterator(__sparse->end()) : iterator(__dense.end())); }
    
    inline const_sparse_iterator sparse_begin() const { return __sparse->begin(); }
    inline       sparse_iterator sparse_begin()       { return __sparse->begin(); }
    inline const_sparse_iterator sparse_end() const { return __sparse->end(); }
    inline       sparse_iterator sparse_end()       { return __sparse->end(); }
    
    inline const_dense_iterator dense_begin() const { return __dense.begin(); }
    inline       dense_iterator dense_begin()       { return __dense.begin(); }
    inline const_dense_iterator dense_end() const { return __dense.end(); }
    inline       dense_iterator dense_end()       { return __dense.end(); }
    
    
    void swap(FeatureVector& x)
    { 
      __dense.swap(x.__dense);
      std::swap(__sparse, x.__sparse);
    }
    
    Tp sum() const
    {
      return (__sparse
	      ? __sum_aux(__sparse->begin(), __sparse->end())
	      : __sum_aux(__dense.begin(), __dense.end()));
    }

    template <typename Iterator>
    Tp __sum_aux(Iterator first, Iterator last) const
    {
      Tp __sum = Tp();
      for (/**/; first != last; ++ first)
	__sum += first->second;
      return __sum;
    }
    
  public:
    // comparison
    friend
    bool operator==(const FeatureVector& x, const FeatureVector& y)
    {
      return (x.__sparse && y.__sparse && *x.__sparse == *y.__sparse) || (x.__sparse == 0 && y.__sparse == 0 && x.__dense == y.__dense);
    }

    friend
    bool operator!=(const FeatureVector& x, const FeatureVector& y)
    {
      return ! (x == y);
    }
    
    friend
    bool operator<(const FeatureVector& x, const FeatureVector& y)
    {
#if 1
      if (x.__sparse && y.__sparse)
	return *x.__sparse < *y.__sparse;
      else if (x.__sparse == 0 && y.__sparse == 0)
	return x.__dense < y.__dense;
      else
	return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
#endif

      //return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
#if 0
      if (x.__sparse)
	return (y.__sparse
		? *x.__sparse < *y.sparse
		: std::lexicographical_compare(x.__sparse->begin(), x.__sparse->end(),
					       y.__dense.begin(), y.__dense.end()));
      else
	return (y.__sparse
		? std::lexicographical_compare(x.__dense.begin(), x.__dense.end(),
					       y.__sparse->begin(), y.__sparse->end())
		: x.__dense < y.__dense);
#endif
    }

    friend
    bool operator<=(const FeatureVector& x, const FeatureVector& y)
    {
      return ! (y < x);
    }

    friend
    bool operator>(const FeatureVector& x, const FeatureVector& y)
    {
      return y < x;
    }

    friend
    bool operator>=(const FeatureVector& x, const FeatureVector& y)
    {
      return ! (x < y);
    }


  public:
    
    template <typename T, typename A>
    friend
    std::ostream& operator<<(std::ostream& os, const FeatureVector<T,A>& x);
    
    template <typename T, typename A>
    friend
    std::istream& operator>>(std::istream& is, FeatureVector<T,A>& x);

  private:
    template <typename O, typename T>
    struct __apply_unary : public O
    {
      __apply_unary(const T& x) : const_value(x) {}
      
      template <typename Value>
      void operator()(Value& value) const
      {
	value.second = O::operator()(value.second, const_value);
      }

      T const_value;
    };
    
  public:
    // operators...
    template <typename T>
    self_type& operator+=(const T& x)
    { 
      if (__sparse)
	std::for_each(__sparse->begin(), __sparse->end(), __apply_unary<std::plus<Tp>, T>(x));
      else
	std::for_each(__dense.begin(), __dense.end(), __apply_unary<std::plus<Tp>, T>(x));
      return *this;
    }

    template <typename T>
    self_type& operator-=(const T& x)
    { 
      if (__sparse)
	std::for_each(__sparse->begin(), __sparse->end(), __apply_unary<std::minus<Tp>, T>(x));
      else
	std::for_each(__dense.begin(), __dense.end(), __apply_unary<std::minus<Tp>, T>(x));
      return *this;
    }
    
    template <typename T>
    self_type& operator*=(const T& x)
    { 
      if (x == T())
	clear();
      else if (__sparse)
	std::for_each(__sparse->begin(), __sparse->end(), __apply_unary<std::multiplies<Tp>, T>(x));
      else
	std::for_each(__dense.begin(), __dense.end(), __apply_unary<std::multiplies<Tp>, T>(x));
      return *this;
    }
    
    template <typename T>
    self_type& operator/=(const T& x)
    {
      if (__sparse)
	std::for_each(__sparse->begin(), __sparse->end(), __apply_unary<std::divides<Tp>, T>(x));
      else
	std::for_each(__dense.begin(), __dense.end(), __apply_unary<std::divides<Tp>, T>(x));
      return *this;
    }
    
    template <typename T, typename A>
    self_type& operator+=(const FeatureVector<T,A>& x)
    {
      typedef FeatureVector<T,A> another_type;
      
      if (x.empty())
	return *this;
      else if (empty()) {
	assign(x);
	return *this;
      }
      
      if (__sparse || x.sparse()) {
	if (! __sparse) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
	
	if (x.sparse())
	  plus_equal_ordered(*__sparse, x.__sparse->begin(), x.__sparse->end());
	else
	  plus_equal_ordered(*__sparse, x.__dense.begin(), x.__dense.end());
      } else {
	dense_vector_type dense_new;

	typename dense_vector_type::const_iterator iter1     = __dense.begin();
	typename dense_vector_type::const_iterator iter1_end = __dense.end();
	
	typename another_type::dense_vector_type::const_iterator iter2     = x.__dense.begin();
	typename another_type::dense_vector_type::const_iterator iter2_end = x.__dense.end();

	while (iter1 != iter1_end && iter2 != iter2_end) {
	  if (iter1->first < iter2->first) {
	    dense_new.insert(dense_new.end(), *iter1);
	    ++ iter1;
	  } else if (iter2->first < iter1->first) {
	    dense_new.insert(dense_new.end(), *iter2);
	    ++ iter2;
	  } else {
	    const Tp value = iter1->second + iter2->second;
	    if (value != Tp())
	      dense_new.insert(dense_new.end(), std::make_pair(iter1->first, value));
	    ++ iter1;
	    ++ iter2;
	  }
	}
	
	if (iter1 != iter1_end)
	  dense_new.insert(iter1, iter1_end);
	if (iter2 != iter2_end)
	  dense_new.insert(iter2, iter2_end);
	
	__dense.swap(dense_new);
	
	if (__dense.size() > __dense_size) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
      }
      
      return *this;
    }

    template <typename T, typename A>
    self_type& operator+=(const FeatureVectorLinear<T,A>& x)
    {
      typedef FeatureVectorLinear<T,A> another_type;
      
      if (x.empty())
	return *this;
      else if (empty()) {
	assign(x);
	return *this;
      }
      
      if (__sparse || x.size() > __dense_size) {
	if (! __sparse) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
	
	plus_equal_ordered(*__sparse, x.begin(), x.end());
      } else {
	dense_vector_type dense_new;
	
	typename dense_vector_type::const_iterator iter1     = __dense.begin();
	typename dense_vector_type::const_iterator iter1_end = __dense.end();
	
	typename another_type::const_iterator iter2     = x.begin();
	typename another_type::const_iterator iter2_end = x.end();

	while (iter1 != iter1_end && iter2 != iter2_end) {
	  if (iter1->first < iter2->first) {
	    dense_new.insert(dense_new.end(), *iter1);
	    ++ iter1;
	  } else if (iter2->first < iter1->first) {
	    dense_new.insert(dense_new.end(), *iter2);
	    ++ iter2;
	  } else {
	    const Tp value = iter1->second + iter2->second;
	    if (value != Tp())
	      dense_new.insert(dense_new.end(), std::make_pair(iter1->first, value));
	    ++ iter1;
	    ++ iter2;
	  }
	}
	
	if (iter1 != iter1_end)
	  dense_new.insert(iter1, iter1_end);
	if (iter2 != iter2_end)
	  dense_new.insert(iter2, iter2_end);
	
	__dense.swap(dense_new);
	
	if (__dense.size() > __dense_size) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
      }
      
      return *this;
    }

    template <typename T, typename A>
    self_type& operator+=(const FeatureVectorUnordered<T,A>& x)
    {
      typedef FeatureVectorUnordered<T,A> another_type;
      
      if (x.empty())
	return *this;
      else if (empty()) {
	assign(x);
	return *this;
      }
      
      if (__sparse || x.size() > __dense_size) {
	if (! __sparse) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
	
	plus_equal_unordered(*__sparse, x.begin(), x.end());
      } else {
	plus_equal_unordered(__dense, x.begin(), x.end());
	
	if (__dense.size() > __dense_size) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
      }
      
      return *this;
    }
    
    template <typename T, typename A>
    self_type& operator-=(const FeatureVector<T,A>& x)
    {
      typedef FeatureVector<T,A> another_type;
      
      if (x.empty()) return *this;
      
      if (__sparse || x.sparse()) {
	if (! __sparse) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}

	if (x.sparse())
	  minus_equal_ordered(*__sparse, x.__sparse->begin(), x.__sparse->end());
	else
	  minus_equal_ordered(*__sparse, x.__dense.begin(), x.__dense.end());
      } else {
	dense_vector_type dense_new;

	typename dense_vector_type::const_iterator iter1     = __dense.begin();
	typename dense_vector_type::const_iterator iter1_end = __dense.end();
	
	typename another_type::dense_vector_type::const_iterator iter2     = x.__dense.begin();
	typename another_type::dense_vector_type::const_iterator iter2_end = x.__dense.end();

	while (iter1 != iter1_end && iter2 != iter2_end) {
	  if (iter1->first < iter2->first) {
	    dense_new.insert(dense_new.end(), *iter1);
	    ++ iter1;
	  } else if (iter2->first < iter1->first) {
	    dense_new.insert(dense_new.end(), std::make_pair(iter2->first, -Tp(iter2->second)));
	    ++ iter2;
	  } else {
	    const Tp value = iter1->second - iter2->second;
	    if (value != Tp())
	      dense_new.insert(dense_new.end(), std::make_pair(iter1->first, value));
	    ++ iter1;
	    ++ iter2;
	  }
	}
	
	if (iter1 != iter1_end)
	  dense_new.insert(iter1, iter1_end);
	for (/**/; iter2 != iter2_end; ++ iter2)
	  dense_new.insert(dense_new.end(), std::make_pair(iter2->first, -Tp(iter2->second)));
	
	__dense.swap(dense_new);
	
	if (__dense.size() > __dense_size) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
      }
      
      return *this;
    }

    template <typename T, typename A>
    self_type& operator-=(const FeatureVectorLinear<T,A>& x)
    {
      typedef FeatureVectorLinear<T,A> another_type;
      
      if (x.empty()) return *this;
      
      if (__sparse || x.size() > __dense_size) {
	if (! __sparse) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
	
	minus_equal_ordered(*__sparse, x.begin(), x.end());
      } else {
	dense_vector_type dense_new;

	typename dense_vector_type::const_iterator iter1     = __dense.begin();
	typename dense_vector_type::const_iterator iter1_end = __dense.end();
	
	typename another_type::const_iterator iter2     = x.begin();
	typename another_type::const_iterator iter2_end = x.end();
	
	while (iter1 != iter1_end && iter2 != iter2_end) {
	  if (iter1->first < iter2->first) {
	    dense_new.insert(dense_new.end(), *iter1);
	    ++ iter1;
	  } else if (iter2->first < iter1->first) {
	    dense_new.insert(dense_new.end(), std::make_pair(iter2->first, -Tp(iter2->second)));
	    ++ iter2;
	  } else {
	    const Tp value = iter1->second - iter2->second;
	    if (value != Tp())
	      dense_new.insert(dense_new.end(), std::make_pair(iter1->first, value));
	    ++ iter1;
	    ++ iter2;
	  }
	}
	
	if (iter1 != iter1_end)
	  dense_new.insert(iter1, iter1_end);
	for (/**/; iter2 != iter2_end; ++ iter2)
	  dense_new.insert(dense_new.end(), std::make_pair(iter2->first, -Tp(iter2->second)));
	
	__dense.swap(dense_new);
	
	if (__dense.size() > __dense_size) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
      }
      
      return *this;
    }

    template <typename T, typename A>
    self_type& operator-=(const FeatureVectorUnordered<T,A>& x)
    {
      typedef FeatureVectorUnordered<T,A> another_type;
      
      if (x.empty()) return *this;
      
      if (__sparse || x.size() > __dense_size) {
	if (! __sparse) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
	
	minus_equal_unordered(*__sparse, x.begin(), x.end());
      } else {
	minus_equal_unordered(__dense, x.begin(), x.end());
	
	if (__dense.size() > __dense_size) {
	  __sparse = new sparse_vector_type(__dense.begin(), __dense.end());
	  __dense.clear();
	}
      }
      
      return *this;
    }
    
    template <typename T, typename A>
    self_type& operator*=(const FeatureVector<T,A>& x)
    {
      typedef FeatureVector<T,A> another_type;
      typedef std::pair<feature_type, data_type> feat_type;
      typedef std::vector<feat_type, std::allocator<feat_type> > feat_set_type;
      
      if (empty() || x.empty()) {
	clear();
	return *this;
      }
      
      feat_set_type feats;
      
      const_iterator iter1     = lower_bound(x.begin()->first);
      const_iterator iter1_end = end();
      
      typename another_type::const_iterator iter2     = (iter1 != iter1_end ? x.lower_bound(iter1->first) : x.begin());
      typename another_type::const_iterator iter2_end = x.end();
      
      while (iter1 != iter1_end && iter2 != iter2_end) {
	if (iter1->first < iter2->first)
	  ++ iter1;
	else if (iter2->first < iter1->first)
	  ++ iter2;
	else {
	  const Tp value = iter1->second * iter2->second;
	  if (value != Tp())
	    feats.push_back(std::make_pair(iter1->first, value));
	  
	  ++ iter1;
	  ++ iter2;
	}
      }
      
      assign(feats.begin(), feats.end());
      
      return *this;
    }
    
    
    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVector<T1,A1> operator+(const FeatureVector<T1,A1>& x, const FeatureVector<T2,A2>& y);

    
    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVector<T1,A1> operator-(const FeatureVector<T1,A1>& x, const FeatureVector<T2,A2>& y);

    
    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVector<T1,A1> operator*(const FeatureVector<T1,A1>& x, const FeatureVector<T2,A2>& y);

  private:
    template <typename Container, typename Iterator1, typename Iterator2>
    void plus_equal_ordered(Container& container, Iterator1 first1, Iterator1 last1, Iterator2 first2, Itertor2 last2)
    {
      
      
    }

    template <typename Container, typename Iterator>
    void plus_equal_ordered(Container& container, Iterator first, Iterator last)
    {
      typename Container::iterator hint = container.begin();
      
      for (/**/; first != last && hint != container.end(); ++ first) {
	std::pair<typename Container::iterator, bool> result = container.insert(*first);
	
	hint = result.first;
	++ hint;
	
	if (! result.second) {
	  result.first->second += first->second;
	  
	  if (result.first->second == Tp()) {
	    container.erase(result.first);
	    hint = container.begin();
	  }
	}
      }
      
      for (/**/; first != last; ++ first)
	container.insert(container.end(), *first);
    }
    
    template <typename Container, typename Iterator>
    void plus_equal_unordered(Container& container, Iterator first, Iterator last)
    {
      for (/**/; first != last; ++ first) {
	std::pair<typename Container::iterator, bool> result = container.insert(*first);
	
	if (! result.second) {
	  result.first->second += first->second;
	  
	  if (result.first->second == Tp())
	    container.erase(result.first);
	}
      }
    }

    template <typename Container, typename Iterator>
    void minus_equal_ordered(Container& container, Iterator first, Iterator last)
    {
      typename Container::iterator hint = container.begin();
      
      for (/**/; first != last && hint != container.end(); ++ first) {
	std::pair<typename Container::iterator, bool> result = container.insert(std::make_pair(first->first, -Tp(first->second)));

	hint = result.first;
	++ hint;
	
	if (! result.second) {
	  result.first->second -= first->second;
	  
	  if (result.first->second == Tp()) {
	    container.erase(result.first);
	    hint = container.begin();
	  }
	}
      }
      
      for (/**/; first != last; ++ first)
	container.insert(container.end(), std::make_pair(first->first, -Tp(first->second)));
    }
    
    template <typename Container, typename Iterator>
    void minus_equal_unordered(Container& container, Iterator first, Iterator last)
    {
      for (/**/; first != last; ++ first) {
	std::pair<typename Container::iterator, bool> result = container.insert(std::make_pair(first->first, -Tp(first->second)));
	
	if (! result.second) {
	  result.first->second -= first->second;
	  
	  if (result.first->second == Tp())
	    container.erase(result.first);
	}
      }
    }

  public:
    dense_vector_type   __dense;
    sparse_vector_type* __sparse;
  };
  
  template <typename T1, typename A1, typename T2>
  inline
  FeatureVector<T1,A1> operator+(const FeatureVector<T1,A1>& x, const T2& y)
  {
    FeatureVector<T1,A1> features(x);
    features += y;
    return features;
  }

  template <typename T2, typename T1, typename A1>
  inline
  FeatureVector<T1,A1> operator+(const T2& x, const FeatureVector<T1,A1>& y)
  {
    FeatureVector<T1,A1> features(y);
    features += x;
    return features;
  }

  template <typename T1, typename A1, typename T2>
  inline
  FeatureVector<T1,A1> operator*(const FeatureVector<T1,A1>& x, const T2& y)
  {
    if (y == T2()) return FeatureVector<T1,A1>();
    
    FeatureVector<T1,A1> features(x);
    features *= y;
    return features;
  }

  template <typename T2, typename T1, typename A1>
  inline
  FeatureVector<T1,A1> operator*(const T2& x, const FeatureVector<T1,A1>& y)
  {
    if (x == T2()) return FeatureVector<T1,A1>();
    
    FeatureVector<T1,A1> features(y);
    features *= x;
    return features;
  }

  template <typename T1, typename A1, typename T2>
  inline
  FeatureVector<T1,A1> operator-(const FeatureVector<T1,A1>& x, const T2& y)
  {
    FeatureVector<T1,A1> features(x);
    features -= y;
    return features;
  }
  
  template <typename T1, typename A1, typename T2>
  inline
  FeatureVector<T1,A1> operator/(const FeatureVector<T1,A1>& x, const T2& y)
  {
    FeatureVector<T1,A1> features(x);
    features /= y;
    return features;
  }

  
  template <typename T1, typename A1, typename T2, typename A2>
  inline
  FeatureVector<T1,A1> operator+(const FeatureVector<T1,A1>& x, const FeatureVector<T2,A2>& y)
  {
    typedef FeatureVector<T1,A1> left_type;
    typedef FeatureVector<T2,A2> right_type;

    if (y.empty())
      return x;
    else if (x.empty())
      return y;
    else {
      left_type features(x);
      
      features += y;
      
      return features;
    }
  }

  template <typename T1, typename A1, typename T2, typename A2>
  inline
  FeatureVector<T1,A1> operator-(const FeatureVector<T1,A1>& x, const FeatureVector<T2,A2>& y)
  {
    typedef FeatureVector<T1,A1> left_type;
    typedef FeatureVector<T2,A2> right_type;

    if (y.empty())
      return x;
    else if (x.empty()) {
      left_type features;
      
      typename right_type::const_iterator iter2_end = y.end();
      for (typename right_type::const_iterator iter2 = y.begin(); iter2 != iter2_end; ++ iter2)
	features.insert(features.end(), std::make_pair(iter2->first, - T1(iter2->second)));
      
      return features;
    } else {
      left_type features(x);
      
      features -= y;
      
      return features;
    }
  }

  template <typename T1, typename A1, typename T2, typename A2>
  inline
  FeatureVector<T1,A1> operator*(const FeatureVector<T1,A1>& x, const FeatureVector<T2,A2>& y)
  {
    typedef FeatureVector<T1,A1> left_type;
    typedef FeatureVector<T2,A2> right_type;
    
    if (x.empty() || y.empty())
      return left_type();
    
    left_type features;
    
    typename left_type::const_iterator iter1     = x.lower_bound(y.begin()->first);
    typename left_type::const_iterator iter1_end = x.end();

    typename right_type::const_iterator iter2     = (iter1 != iter1_end ? y.lower_bound(iter1->first) : y.begin());
    typename right_type::const_iterator iter2_end = y.end();
    
    while (iter1 != iter1_end && iter2 != iter2_end) {
      if (iter1->first < iter2->first)
	++ iter1;
      else if (iter2->first < iter1->first)
	++ iter2;
      else {
	const T1 value = iter1->second * iter2->second;
	if (value != T1())
	  features.insert(features.end(), std::make_pair(iter1->first, value));
	
	++ iter1;
	++ iter2;
      }
    }
    
    return features;
  }


  template <typename T, typename A>
  inline
  std::ostream& operator<<(std::ostream& os, const FeatureVector<T,A>& x)
  {
    typename FeatureVector<T,A>::const_iterator iter_end = x.end();
    for (typename FeatureVector<T,A>::const_iterator iter = x.begin(); iter != iter_end; ++ iter)
      if (! iter->first.empty() && iter->second != T())
	os << iter->first << ' ' << iter->second << '\n';
    
    return os;
  }

  template <typename T, typename A>
  inline
  std::istream& operator>>(std::istream& is, FeatureVector<T,A>& x)
  {
    x.clear();
    
    std::string feature;
    T value;
    while ((is >> feature) && (is >> value))
      if (value != T())
	x[feature] = value;
    
    return is;
  }
  

};

namespace std
{
  template <typename T, typename A>
  inline
  void swap(cicada::FeatureVector<T,A>& x, cicada::FeatureVector<T, A>& y)
  {
    x.swap(y);
  }
};

#include <cicada/weight_vector.hpp>
#include <cicada/feature_vector_compact.hpp>
#include <cicada/feature_vector_linear.hpp>
#include <cicada/feature_vector_unordered.hpp>

namespace cicada
{
  template <typename T, typename A>
  inline
  void FeatureVector<T,A>::assign(const FeatureVectorCompact& x)
  {
    typedef std::pair<feature_type, data_type> feat_type;
    typedef std::vector<feat_type, std::allocator<feat_type> > feat_set_type;
    
    feat_set_type feats(x.begin(), x.end());
    
    assign(feats.begin(), feats.end());
  }
};  

#endif

