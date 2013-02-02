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

#include <utils/compact_map.hpp>
#include <utils/hashmurmur.hpp>

namespace cicada
{
  
  // forward declaration...
  template <typename Tp, typename Alloc >
  class WeightVector;

  class FeatureVectorCompact;

  template <typename Tp, typename Alloc >
  class FeatureVectorLinear;

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

    typedef typename Alloc::template rebind<value_type>::other alloc_type;
    typedef typename utils::compact_map<key_type, data_type,
					utils::unassigned<key_type>, utils::deleted<key_type>,
					utils::hashmurmur<size_t>, std::equal_to<key_type>,
					alloc_type> vector_type;
    
    typedef FeatureVector<Tp, Alloc> self_type;
    
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef typename vector_type::const_iterator const_iterator;
    typedef typename vector_type::iterator             iterator;
    
    typedef typename vector_type::const_reference const_reference;
    typedef typename vector_type::reference       reference;
    typedef typename vector_type::pointer         pointer;
    
  public:
    FeatureVector(size_type hint=8) : __vector(hint) { rehash(hint); }
    FeatureVector(const FeatureVector<Tp,Alloc>& x) : __vector(x.__vector) {}
    template <typename T, typename A>
    FeatureVector(const FeatureVector<T,A>& x) : __vector(x.size()) { assign(x); }
    template <typename Iterator>
    FeatureVector(Iterator first, Iterator last) : __vector() { assign(first, last); }
    FeatureVector(const FeatureVectorCompact& x) : __vector() { assign(x); }
    template <typename T, typename A>
    FeatureVector(const FeatureVectorLinear<T,A>& x): __vector(x.size()) { assign(x); } 
    
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
    
  public:
    void assign(const FeatureVector<Tp,Alloc>& x)
    {
      __vector = x.__vector;
    }
    
    template <typename T, typename A>
    void assign(const FeatureVector<T,A>& x)
    {
      __vector.clear();
      __vector.rehash(x.size());
      __vector.insert(x.begin(), x.end());
    }

    void assign(const FeatureVectorCompact& x);

    template <typename T, typename A>
    void assign(const FeatureVectorLinear<T,A>& x)
    {
      __vector.clear();
      __vector.rehash(x.size());
      __vector.insert(x.begin(), x.end());
    }

    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      __vector.clear();
      __vector.insert(first, last);
    }
    
    template <typename Iterator>
    void insert(Iterator first, Iterator last)
    {
      __vector.insert(first, last);
    }
    
    void insert(iterator iter, const value_type& x)
    {
      __vector.insert(iter, x);
    }
    
    void insert(const value_type& x)
    {
      __vector.insert(x);
    }

    template <typename T, typename A>
    FeatureVector& intersect(const FeatureVector<T,A>& x)
    {
      if (empty()) 
	return *this;
      
      if (x.empty())
	clear();
      else {
	vector_type vector_new;
	vector_new.rehash(utils::bithack::max(__vector.size(), x.size()));
	
	intersect(vector_new, __vector, x.begin(), x.end());
	
	__vector.swap(vector_new);
      }
      
      return *this; 
    }

    FeatureVector& intersect(const FeatureVectorCompact& x);

    template <typename T, typename A>
    FeatureVector& intersect(const FeatureVectorLinear<T,A>& x)
    {
      if (empty()) 
	return *this;
      
      if (x.empty())
	clear();
      else {
	vector_type vector_new;
	vector_new.rehash(utils::bithack::max(__vector.size(), x.size()));
	
	intersect(vector_new, __vector, x.begin(), x.end());
	
	__vector.swap(vector_new);
      }
      
      return *this; 
    }

    template <typename T, typename A>
    FeatureVector& intersect_absmax(const FeatureVector<T,A>& x)
    {
      if (empty()) 
	operator=(x);
      else if (! x.empty()) {
	vector_type vector_new;
	vector_new.rehash(utils::bithack::max(__vector.size(), x.size()));
	
	intersect_absmax(vector_new, __vector, x.begin(), x.end());
	
	__vector.swap(vector_new);
      }
      
      return *this; 
    }

    FeatureVector& intersect_absmax(const FeatureVectorCompact& x);

    template <typename T, typename A>
    FeatureVector& intersect_absmax(const FeatureVectorLinear<T,A>& x)
    {
      if (empty()) 
	operator=(x);
      else if (! x.empty()) {
	vector_type vector_new;
	vector_new.rehash(utils::bithack::max(__vector.size(), x.size()));
	
	intersect_absmax(vector_new, __vector, x.begin(), x.end());
	
	__vector.swap(vector_new);
      }
      
      return *this; 
    }
    
    template <typename T, typename A>
    FeatureVector& intersect_absmin(const FeatureVector<T,A>& x)
    {
      if (empty()) 
	return *this;
      
      if (x.empty())
	clear();
      else {
	vector_type vector_new;
	vector_new.rehash(utils::bithack::max(__vector.size(), x.size()));
	
	intersect_absmin(vector_new, __vector, x.begin(), x.end());
	
	__vector.swap(vector_new);
      }
      
      return *this; 
    }

    FeatureVector& intersect_absmin(const FeatureVectorCompact& x);

    template <typename T, typename A>
    FeatureVector& intersect_absmin(const FeatureVectorLinear<T,A>& x)
    {
      if (empty()) 
	return *this;
      
      if (x.empty())
	clear();
      else {
	vector_type vector_new;
	vector_new.rehash(utils::bithack::max(__vector.size(), x.size()));
	
	intersect_absmin(vector_new, __vector, x.begin(), x.end());
	
	__vector.swap(vector_new);
      }
      
      return *this; 
    }

    template <typename T, typename A, typename Prefix>
    void update(const FeatureVector<T,A>& x, const Prefix& prefix)
    {
      if (empty())
	assign(x);
      else {
	erase_prefix(prefix);
	operator+=(x);
      }
    }
    
    template <typename Prefix>
    void update(const FeatureVectorCompact& x, const Prefix& prefix)
    {
      if (empty())
	assign(x);
      else {
	erase_prefix(prefix);
	operator+=(x);
      }
    }
    
    template <typename T, typename A, typename Prefix>
    void update(const FeatureVectorLinear<T,A>& x, const Prefix& prefix)
    {
      if (empty())
	assign(x);
      else {
	erase_prefix(prefix);
	operator+=(x);
      }
    }
    
  private:
    template <typename Container, typename Original, typename Iterator>
    static inline
    void intersect(Container& container, const Original& orig, Iterator first, Iterator last)
    {
      typedef typename Original::value_type::second_type                       value1_type;
      typedef typename std::iterator_traits<Iterator>::value_type::second_type value2_type;
      
      for (/**/; first != last; ++ first) {
	typename Original::const_iterator iter = orig.find(first->first);
	
	if (iter == orig.end()) continue;
	
	if (iter->second > value1_type() && first->second > value2_type()) {
	  const value1_type value(std::min(iter->second, static_cast<value1_type>(first->second)));
	  
	  container.insert(std::make_pair(iter->first, value));
	} else if (iter->second < value1_type() && first->second < value2_type()) {
	  const value1_type value(std::max(iter->second, static_cast<value1_type>(first->second)));
	  
	  container.insert(std::make_pair(iter->first, value));
	}
      }
    }

    template <typename Container, typename Original, typename Iterator>
    static inline
    void intersect_absmax(Container& container, const Original& orig, Iterator first, Iterator last)
    {
      typedef typename Original::value_type::second_type                       value1_type;
      typedef typename std::iterator_traits<Iterator>::value_type::second_type value2_type;
      
      for (/**/; first != last; ++ first) {
	typename Original::const_iterator iter = orig.find(first->first);
	
	if (iter == orig.end()) {
	  if (first->second != value2_type())
	    container.insert(*first);
	} else {
	  const value1_type abs1 = (iter->second < value1_type()  ? - iter->second : iter->second);
	  const value2_type abs2 = (first->second < value2_type() ? - first->second : first->second);
	  
	  if (abs1 > abs2) {
	    if (iter->second != value1_type())
	      container.insert(*iter);
	  } else {
	    if (first->second != value2_type())
	      container.insert(*first);
	  }
	}
      }
    }
    
    template <typename Container, typename Original, typename Iterator>
    static inline
    void intersect_absmin(Container& container, const Original& orig, Iterator first, Iterator last)
    {
      typedef typename Original::value_type::second_type                       value1_type;
      typedef typename std::iterator_traits<Iterator>::value_type::second_type value2_type;
      
      for (/**/; first != last; ++ first) {
	typename Original::const_iterator iter = orig.find(first->first);
	
	if (iter == orig.end()) continue;
	if (iter->second == value1_type() || first->second == value2_type()) continue;
	
	const value1_type abs1 = (iter->second < value1_type()  ? - iter->second : iter->second);
	const value2_type abs2 = (first->second < value2_type() ? - first->second : first->second);
	
	if (abs1 < abs2)
	  container.insert(*iter);
	else
	  container.insert(*first);
      }
    }
    
  public:

    size_type size() const { return __vector.size(); }
    bool empty() const { return __vector.empty(); }

    void reserve(size_type x) { __vector.rehash(x); }
    void rehash(size_type x) { __vector.rehash(x); }
    
    void clear()
    {
      __vector.clear();
    }
    
    Tp operator[](const key_type& x) const
    {
      const_iterator iter = __vector.find(x);
      return (iter == __vector.end() ? Tp() : iter->second);
    }
    
    Tp& operator[](const key_type& x)
    {
      return __vector[x];
    }
    
    const_iterator find(const key_type& x) const
    {
      return __vector.find(x);
    }
    
    iterator find(const key_type& x)
    {
      return __vector.find(x);
    }
    
    void erase(const key_type& x)
    {
      __vector.erase(x);
    }

    void erase(iterator x)
    {
      __vector.erase(x.diter);
    }
    
    template <typename Prefix>
    void erase_prefix(const Prefix& prefix)
    {
      for (iterator fiter = __vector.begin(); fiter != __vector.end(); /**/)
	if (fiter->first.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), fiter->first.begin()))
	  __vector.erase(fiter ++);
	else
	  ++ fiter;
    }
    
    inline const_iterator begin() const { return __vector.begin(); }
    inline       iterator begin()       { return __vector.begin(); }
    
    inline const_iterator end() const { return __vector.end(); }
    inline       iterator end()       { return __vector.end(); }
    
    
    void swap(FeatureVector& x)
    { 
      __vector.swap(x.__vector);
    }
    
    Tp sum() const
    {
      return __sum_aux(__vector.begin(), __vector.end());
    }

  private:
    template <typename Iterator>
    Tp __sum_aux(Iterator first, Iterator last) const
    {
      Tp __sum = Tp();
      for (/**/; first != last; ++ first)
	__sum += first->second;
      return __sum;
    }

  private:
    struct __equal_to
    {
      bool operator()(const vector_type& x, const vector_type& y) const
      {
	if (x.size() != y.size()) return false;
	
	typename vector_type::const_iterator iter_end = y.end();
	for (typename vector_type::const_iterator iter = y.begin(); iter != iter_end; ++ iter) {
	  typename vector_type::const_iterator fiter = x.find(iter->first);
	  
	  if (fiter == x.end() || *fiter != *iter)
	    return false;
	}
	
	return true;
      }
    };
    
  public:
    // comparison
    friend
    bool operator==(const FeatureVector& x, const FeatureVector& y)
    {
      return __equal_to()(x.__vector, y.__vector);
    }
    
    friend
    bool operator!=(const FeatureVector& x, const FeatureVector& y)
    {
      return ! (x == y);
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
      std::for_each(__vector.begin(), __vector.end(), __apply_unary<std::plus<Tp>, T>(x));
      return *this;
    }

    template <typename T>
    self_type& operator-=(const T& x)
    { 
      std::for_each(__vector.begin(), __vector.end(), __apply_unary<std::minus<Tp>, T>(x));
      return *this;
    }
    
    template <typename T>
    self_type& operator*=(const T& x)
    { 
      if (x == T())
	clear();
      else
	std::for_each(__vector.begin(), __vector.end(), __apply_unary<std::multiplies<Tp>, T>(x));
      return *this;
    }
    
    template <typename T>
    self_type& operator/=(const T& x)
    {
      std::for_each(__vector.begin(), __vector.end(), __apply_unary<std::divides<Tp>, T>(x));
      return *this;
    }
    
    template <typename T, typename A>
    self_type& operator+=(const FeatureVector<T,A>& x)
    {
      if (x.empty())
	return *this;
      else if (empty()) {
	assign(x);
	return *this;
      } else {
	__vector.rehash(__vector.occupied_count() + x.size());
	plus_equal(__vector, x.begin(), x.end());
	return *this;
      }
    }

    self_type& operator+=(const FeatureVectorCompact& x);

    template <typename T, typename A>
    self_type& operator+=(const FeatureVectorLinear<T,A>& x)
    {
      if (x.empty())
	return *this;
      else if (empty()) {
	assign(x);
	return *this;
      } else {
	__vector.rehash(__vector.occupied_count() + x.size());
	plus_equal(__vector, x.begin(), x.end());
	return *this;
      }
    }

    template <typename T, typename A>
    self_type& operator-=(const FeatureVector<T,A>& x)
    {
      if (x.empty()) return *this;
      
      __vector.rehash(__vector.occupied_count() + x.size());
      minus_equal(__vector, x.begin(), x.end());
      
      return *this;
    }
    
    self_type& operator-=(const FeatureVectorCompact& x);

    template <typename T, typename A>
    self_type& operator-=(const FeatureVectorLinear<T,A>& x)
    {
      if (x.empty()) return *this;
      
      __vector.rehash(__vector.occupied_count() + x.size());
      minus_equal(__vector, x.begin(), x.end());
      
      return *this;
    }

    template <typename T, typename A>
    self_type& operator*=(const FeatureVector<T,A>& x)
    {
      if (empty() || x.empty()) {
	clear();
	return *this;
      } else {
	vector_type vector_new;
	vector_new.rehash(utils::bithack::max(__vector.size(), x.size()));
	
	multiply_equal(vector_new, __vector, x.begin(), x.end());
	
	__vector.swap(vector_new);
	
	return *this;
      }
    }

    self_type& operator*=(const FeatureVectorCompact& x);

    template <typename T, typename A>
    self_type& operator*=(const FeatureVectorLinear<T,A>& x)
    {
      if (empty() || x.empty()) {
	clear();
	return *this;
      } else {
	vector_type vector_new;
	vector_new.rehash(utils::bithack::max(__vector.size(), x.size()));
	
	multiply_equal(vector_new, __vector, x.begin(), x.end());
	
	__vector.swap(vector_new);
	
	return *this;
      }
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

    template <typename T1, typename A1>
    friend
    FeatureVector<T1,A1> operator+(const FeatureVector<T1,A1>& x, const FeatureVectorCompact& y);
    
    template <typename T1, typename A1>
    friend
    FeatureVector<T1,A1> operator-(const FeatureVector<T1,A1>& x, const FeatureVectorCompact& y);
    
    template <typename T1, typename A1>
    friend
    FeatureVector<T1,A1> operator*(const FeatureVector<T1,A1>& x, const FeatureVectorCompact& y);

    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVector<T1,A1> operator+(const FeatureVector<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y);
    
    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVector<T1,A1> operator-(const FeatureVector<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y);
    
    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVector<T1,A1> operator*(const FeatureVector<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y);

  private:
    template <typename Container, typename Iterator>
    static inline
    void plus_equal(Container& container, Iterator first, Iterator last)
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
    static inline
    void minus_equal(Container& container, Iterator first, Iterator last)
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
    
    template <typename Container, typename Original, typename Iterator>
    static inline
    void multiply_equal(Container& container, const Original& orig, Iterator first, Iterator last)
    {
      for (/**/; first != last; ++ first) {
	typename Original::const_iterator iter = orig.find(first->first);
	
	if (iter == orig.end()) continue;
	
	const Tp value(iter->second * first->second);
	
	if (value != Tp())
	  container.insert(std::make_pair(first->first, value));
      }
    }

  public:
    vector_type __vector;
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
      left_type features(y.size());
      
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
    
    left_type features(utils::bithack::max(x.size(), y.size()));

    left_type::multiply_equal(features, x, y.begin(), y.end());
    
    return features;
  }

  template <typename T1, typename A1, typename T2, typename A2>
  inline
  FeatureVector<T1,A1> operator+(const FeatureVector<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y)
  {
    typedef FeatureVector<T1,A1> left_type;
    typedef FeatureVectorLinear<T2,A2> right_type;

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
  FeatureVector<T1,A1> operator-(const FeatureVector<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y)
  {
    typedef FeatureVector<T1,A1> left_type;
    typedef FeatureVectorLinear<T2,A2> right_type;

    if (y.empty())
      return x;
    else if (x.empty()) {
      left_type features(y.size());
      
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
  FeatureVector<T1,A1> operator*(const FeatureVector<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y)
  {
    typedef FeatureVector<T1,A1> left_type;
    typedef FeatureVectorLinear<T2,A2> right_type;
    
    if (x.empty() || y.empty())
      return left_type();
    
    left_type features(utils::bithack::max(x.size(), y.size()));

    left_type::multiply_equal(features, x, y.begin(), y.end());
    
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
  
  template <typename T, typename A>
  inline
  FeatureVector<T,A>& FeatureVector<T,A>::intersect(const FeatureVectorCompact& x)
  {
    if (empty()) 
      return *this;
    
    if (x.empty())
      clear();
    else {
      vector_type vector_new;
      vector_new.rehash(__vector.size());
      
      intersect(vector_new, __vector, x.begin(), x.end());
      
      __vector.swap(vector_new);
    }
    
    return *this; 
  }

  template <typename T, typename A>
  inline
  FeatureVector<T,A>& FeatureVector<T,A>::intersect_absmax(const FeatureVectorCompact& x)
  {
    if (empty())
      operator=(x);
    else if (! x.empty()) {
      vector_type vector_new;
      vector_new.rehash(__vector.size());
      
      intersect_absmax(vector_new, __vector, x.begin(), x.end());
      
      __vector.swap(vector_new);
    }
    
    return *this; 
  }
  
  template <typename T, typename A>
  inline
  FeatureVector<T,A>& FeatureVector<T,A>::intersect_absmin(const FeatureVectorCompact& x)
  {
    if (empty())
      return *this;
    
    if (x.empty())
      clear();
    else {
      vector_type vector_new;
      vector_new.rehash(__vector.size());
      
      intersect_absmin(vector_new, __vector, x.begin(), x.end());
      
      __vector.swap(vector_new);
    }
    
    return *this; 
  }


  template <typename T, typename A>
  inline
  FeatureVector<T,A>& FeatureVector<T,A>::operator+=(const FeatureVectorCompact& x)
  {
    if (x.empty())
      return *this;
    else if (empty()) {
      assign(x);
      return *this;
    } else {
      plus_equal(__vector, x.begin(), x.end());
      return *this;
    }
  }

  template <typename T, typename A>
  inline
  FeatureVector<T,A>& FeatureVector<T,A>::operator-=(const FeatureVectorCompact& x)
  {
    if (x.empty()) return *this;

    minus_equal(__vector, x.begin(), x.end());
    
    return *this;
  }

  template <typename T, typename A>
  inline
  FeatureVector<T,A>& FeatureVector<T,A>::operator*=(const FeatureVectorCompact& x)
  {
    if (empty() || x.empty()) {
      clear();
      return *this;
    } else {
      vector_type vector_new;
      vector_new.rehash(__vector.size());
      
      multiply_equal(vector_new, __vector, x.begin(), x.end());
      
      __vector.swap(vector_new);
      
      return *this;
    }
  }
  
  template <typename T1, typename A1>
  inline
  FeatureVector<T1,A1> operator+(const FeatureVector<T1,A1>& x, const FeatureVectorCompact& y)
  {
    typedef FeatureVector<T1,A1> left_type;
    typedef FeatureVectorCompact right_type;

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

  template <typename T1, typename A1>
  inline
  FeatureVector<T1,A1> operator-(const FeatureVector<T1,A1>& x, const FeatureVectorCompact& y)
  {
    typedef FeatureVector<T1,A1> left_type;
    typedef FeatureVectorCompact right_type;

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

  template <typename T1, typename A1>
  inline
  FeatureVector<T1,A1> operator*(const FeatureVector<T1,A1>& x, const FeatureVectorCompact& y)
  {
    typedef FeatureVector<T1,A1> left_type;
    typedef FeatureVectorCompact right_type;
    
    if (x.empty() || y.empty())
      return left_type();
    
    left_type features(x.size());

    left_type::multiply_equal(features, x, y.begin(), y.end());
    
    return features;
  }
};  

#endif

