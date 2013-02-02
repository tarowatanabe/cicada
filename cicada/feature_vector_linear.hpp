// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE_VECTOR_LINEAR__HPP__
#define __CICADA__FEATURE_VECTOR_LINEAR__HPP__ 1

#include <map>
#include <memory>
#include <utility>
#include <algorithm>
#include <iterator>
#include <iostream>

#include <cicada/feature.hpp>

#include <utils/vector_map.hpp>
#include <utils/hashmurmur.hpp>

//
// feature vector, for use as "temporary" linear vector for faster access, w/o sorting
//

namespace cicada
{
  
  // forward declaration...
  template <typename Tp, typename Alloc >
  class FeatureVector;

  class FeatureVectorCompact;
    
  template <typename Tp, typename Alloc >
  class WeightVector;
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class FeatureVectorLinear
  {
  public:
    typedef cicada::Feature feature_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

  private:
    typedef std::pair<feature_type, Tp> __value_type;
    typedef typename Alloc::template rebind<__value_type>::other alloc_type;
    typedef utils::vector_map<feature_type, Tp, std::less<feature_type>,  alloc_type> map_type;
    
  public:
    typedef typename map_type::key_type    key_type;
    typedef typename map_type::data_type   data_type;
    typedef typename map_type::mapped_type mapped_type;
    typedef typename map_type::value_type  value_type;
    
    typedef typename map_type::const_iterator  const_iterator;
    typedef typename map_type::iterator        iterator;
     
    typedef typename map_type::const_reference  const_reference;
    typedef typename map_type::reference        reference;
     
  public:
    typedef FeatureVectorLinear<Tp, Alloc> self_type;
     
  public:
    FeatureVectorLinear() : __map() { }
    FeatureVectorLinear(const self_type& x) : __map(x.__map) { }
    template <typename T, typename A>
    FeatureVectorLinear(const FeatureVector<T,A>& x) : __map() { assign(x); }
    FeatureVectorLinear(const FeatureVectorCompact& x) : __map() { assign(x); }
    template <typename T, typename A>
    FeatureVectorLinear(const FeatureVectorLinear<T,A>& x) : __map() { assign(x); }
    template <typename Iterator>
    FeatureVectorLinear(Iterator first, Iterator last) : __map() { __map.insert(first, last); } 
    
    FeatureVectorLinear& operator=(const self_type& x)
    {
      __map = x.__map;
      return *this;
    }
    
    template <typename T, typename A>
    FeatureVectorLinear& operator=(const FeatureVector<T,A>& x)
    {
      assign(x);
      return *this;
    }
     
    FeatureVectorLinear& operator=(const FeatureVectorCompact& x)
    {
      assign(x);
      return *this;
    }

    template <typename T, typename A>
    FeatureVectorLinear& operator=(const FeatureVectorLinear<T,A>& x)
    {
      assign(x);
      return *this;
    }

  public:
    size_type size() const { return __map.size(); }
    bool empty() const { return __map.empty(); }
     
    void assign(const self_type& x)
    {
      __map = x.__map;
    }

    template <typename T, typename A>
    void assign(const FeatureVector<T,A>& x)
    {
      __map.clear();
      __map.insert(x.begin(), x.end());
    }

    void assign(const FeatureVectorCompact& x);

    template <typename T, typename A>
    void assign(const FeatureVectorLinear<T,A>& x)
    {
      __map.clear();
      __map.insert(x.begin(), x.end());
    }
    
    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      __map.clear();
      __map.insert(first, last);
    }
     
    Tp operator[](const key_type& x) const
    {
      const_iterator iter = find(x);
      return (iter == end() ? Tp() : iter->second);
    }
     
    Tp& operator[](const key_type& x)
    {
      return __map[x];
    }

    inline const_iterator begin() const { return __map.begin(); }
    inline       iterator begin()       { return __map.begin(); }
     
    inline const_iterator end() const { return __map.end(); }
    inline       iterator end()       { return __map.end(); }

    reference front() { return __map.front(); }
    const_reference front() const { return __map.front(); }
    reference back() { return __map.back(); }
    const_reference back() const { return __map.back(); }
    
    inline const_iterator find(const key_type& x) const { return __map.find(x); }
    inline       iterator find(const key_type& x)       { return __map.find(x); }
    
    inline const_iterator lower_bound(const key_type& x) const { return __map.lower_bound(x); }
    inline       iterator lower_bound(const key_type& x)       { return __map.lower_bound(x); }
    
    inline const_iterator upper_bound(const key_type& x) const { return __map.upper_bound(x); }
    inline       iterator upper_bound(const key_type& x)       { return __map.upper_bound(x); }

    void erase(const key_type& x) { __map.erase(x); }
     
    void swap(FeatureVectorLinear& x) { __map.swap(x.__map); }
     
    void clear() { __map.clear(); }

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
    //operators...
     
    template <typename T>
    FeatureVectorLinear& operator+=(const T& x)
    {
      std::for_each(begin(), end(), __apply_unary<std::plus<Tp>, T>(x));
      return *this;
    }

    template <typename T>
    FeatureVectorLinear& operator-=(const T& x)
    {
      std::for_each(begin(), end(), __apply_unary<std::minus<Tp>, T>(x));
      return *this;
    }
     
    template <typename T>
    FeatureVectorLinear& operator*=(const T& x)
    {
      if (x == T())
	clear();
      else
	std::for_each(begin(), end(), __apply_unary<std::multiplies<Tp>, T>(x));
      return *this;
    }
    
    template <typename T>
    FeatureVectorLinear& operator/=(const T& x)
    {
      std::for_each(begin(), end(), __apply_unary<std::divides<Tp>, T>(x));
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
	plus_equal(__map, x.begin(), x.end());
	return *this;
      }
    }

    template <typename T, typename A>
    self_type& operator-=(const FeatureVector<T,A>& x)
    {
      minus_equal(__map, x.begin(), x.end());
      return *this;
    }


    template <typename T, typename A>
    self_type& operator*=(const FeatureVector<T,A>& x)
    {
      if (empty() || x.empty()) {
	clear();
	return *this;
      } else {
	map_type map_new;
	multiply_equal(map_new, __map, x.begin(), x.end());
	__map.swap(map_new);
	return *this;
      }
    }


    template <typename T, typename A>
    self_type& operator+=(const FeatureVectorLinear<T,A>& x)
    {
      if (x.empty())
	return *this;
      else if (empty()) {
	assign(x);
	return *this;
      } else {
	plus_equal_ordered(__map, x.begin(), x.end());
	return *this;
      }
    }
    
    template <typename T, typename A>
    self_type& operator-=(const FeatureVectorLinear<T,A>& x)
    {
      minus_equal_ordered(__map, x.begin(), x.end());
      return *this;
    }
    
    template <typename T, typename A>
    self_type& operator*=(const FeatureVectorLinear<T,A>& x)
    {
      if (empty() || x.empty()) {
	clear();
	return *this;
      } else {
	map_type map_new;
	multiply_equal_ordered(map_new, __map.begin(), __map.end(), x.begin(), x.end());
	__map.swap(map_new);
	return *this;
      }
    }
    self_type& operator+=(const FeatureVectorCompact& x);
    self_type& operator-=(const FeatureVectorCompact& x);
    self_type& operator*=(const FeatureVectorCompact& x);

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
    void plus_equal_ordered(Container& container, Iterator first, Iterator last)
    {
      typename Container::iterator hint = container.begin();

      for (/**/; first != last && hint != container.end(); ++ first) {
	std::pair<typename Container::iterator, bool> result = container.insert(*first);
	
	hint = result.first;
        ++ hint;
	
	if (! result.second) {
	  result.first->second += first->second;
	  
	  if (result.first->second == Tp())
	    container.erase(result.first);
	}
      }
      
      if (first != last)
	container.insert(first, last);
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

    template <typename Container, typename Iterator>
    static inline
    void minus_equal_ordered(Container& container, Iterator first, Iterator last)
    {
      typename Container::iterator hint = container.begin();

      for (/**/; first != last && hint != container.end(); ++ first) {
	std::pair<typename Container::iterator, bool> result = container.insert(std::make_pair(first->first, -Tp(first->second)));
	
	hint = result.first;
        ++ hint;
	
	if (! result.second) {
	  result.first->second -= first->second;
	  
	  if (result.first->second == Tp())
	    container.erase(result.first);
	}
      }
      
      for (/**/; first != last; ++ first)
        container.insert(container.end(), std::make_pair(first->first, -Tp(first->second)));
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

    template <typename Container, typename Iterator1, typename Iterator2>
    static inline
    void multiply_equal_ordered(Container& container,
				Iterator1 first1, Iterator1 last1,
				Iterator2 first2, Iterator2 last2)
    {
      while (first1 != last1 && first2 != last2) {
	if (first1->first < first2->first)
	  ++ first1;
	else if (first2->first < first1->first)
	  ++ first2;
	else {
	  const Tp value = first1->second * first2->second;
	  
	  if (value != Tp())
	    container.insert(container.end(), std::make_pair(first1->first, value));
	  
	  ++ first1;
	  ++ first2;
	}
      }
    }
    
    
  public:
    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVectorLinear<T1,A1> operator+(const FeatureVectorLinear<T1,A1>& x, const FeatureVector<T2,A2>& y);
    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVectorLinear<T1,A1> operator-(const FeatureVectorLinear<T1,A1>& x, const FeatureVector<T2,A2>& y);
    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVectorLinear<T1,A1> operator*(const FeatureVectorLinear<T1,A1>& x, const FeatureVector<T2,A2>& y);

    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVectorLinear<T1,A1> operator+(const FeatureVectorLinear<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y);
    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVectorLinear<T1,A1> operator-(const FeatureVectorLinear<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y);
    template <typename T1, typename A1, typename T2, typename A2>
    friend
    FeatureVectorLinear<T1,A1> operator*(const FeatureVectorLinear<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y);
    
    template <typename T1, typename A1>
    friend
    FeatureVectorLinear<T1,A1> operator+(const FeatureVectorLinear<T1,A1>& x, const FeatureVectorCompact& y);
    
    template <typename T1, typename A1>
    friend
    FeatureVectorLinear<T1,A1> operator-(const FeatureVectorLinear<T1,A1>& x, const FeatureVectorCompact& y);
    
    template <typename T1, typename A1>
    friend
    FeatureVectorLinear<T1,A1> operator*(const FeatureVectorLinear<T1,A1>& x, const FeatureVectorCompact& y);
    
    friend
    size_t hash_value(FeatureVectorLinear const& x) { return utils::hashmurmur<size_t>()(x.__map.begin(), x.__map.end(), 0); }
    
    friend bool operator==(const FeatureVectorLinear& x, const FeatureVectorLinear& y) { return x.__map == y.__map; }
    friend bool operator!=(const FeatureVectorLinear& x, const FeatureVectorLinear& y) { return x.__map != y.__map; }
    friend bool operator<(const FeatureVectorLinear& x, const FeatureVectorLinear& y) { return x.__map < y.__map; }
    friend bool operator>(const FeatureVectorLinear& x, const FeatureVectorLinear& y) { return x.__map > y.__map; }
    friend bool operator<=(const FeatureVectorLinear& x, const FeatureVectorLinear& y) { return x.__map <= y.__map; }
    friend bool operator>=(const FeatureVectorLinear& x, const FeatureVectorLinear& y) { return x.__map >= y.__map; }
     
  private:
    map_type __map;
  };
  

  template <typename T1, typename A1, typename T2>
  inline
  FeatureVectorLinear<T1,A1> operator+(const FeatureVectorLinear<T1,A1>& x, const T2& y)
  {
    FeatureVectorLinear<T1,A1> features(x);
    features += y;
    return features;
  }

  template <typename T2, typename T1, typename A1>
  inline
  FeatureVectorLinear<T1,A1> operator+(const T2& x, const FeatureVectorLinear<T1,A1>& y)
  {
    FeatureVectorLinear<T1,A1> features(y);
    features += x;
    return features;
  }

  template <typename T1, typename A1, typename T2>
  inline
  FeatureVectorLinear<T1,A1> operator*(const FeatureVectorLinear<T1,A1>& x, const T2& y)
  {
    if (y == T2()) return FeatureVectorLinear<T1,A1>();
    
    FeatureVectorLinear<T1,A1> features(x);
    features *= y;
    return features;
  }

  template <typename T2, typename T1, typename A1>
  inline
  FeatureVectorLinear<T1,A1> operator*(const T2& x, const FeatureVectorLinear<T1,A1>& y)
  {
    if (x == T2()) return FeatureVectorLinear<T1,A1>();
    
    FeatureVectorLinear<T1,A1> features(y);
    features *= x;
    return features;
  }

  template <typename T1, typename A1, typename T2>
  inline
  FeatureVectorLinear<T1,A1> operator-(const FeatureVectorLinear<T1,A1>& x, const T2& y)
  {
    FeatureVectorLinear<T1,A1> features(x);
    features -= y;
    return features;
  }
  
  template <typename T1, typename A1, typename T2>
  inline
  FeatureVectorLinear<T1,A1> operator/(const FeatureVectorLinear<T1,A1>& x, const T2& y)
  {
    FeatureVectorLinear<T1,A1> features(x);
    features /= y;
    return features;
  }

  template <typename T1, typename A1, typename T2, typename A2>
  inline
  FeatureVectorLinear<T1,A1> operator+(FeatureVectorLinear<T1,A1>& x, const FeatureVector<T2,A2>& y)
  {
    typedef FeatureVectorLinear<T1,A1> self_type;
    
    self_type x_new(x);
    x_new += y;
    return x_new;
  }
  
  template <typename T1, typename A1, typename T2, typename A2>
  inline
  FeatureVectorLinear<T1,A1> operator-(FeatureVectorLinear<T1,A1>& x, const FeatureVector<T2,A2>& y)
  {
    typedef FeatureVectorLinear<T1,A1> self_type;
    
    self_type x_new(x);
    x_new -= y;
    return x_new;    
  }
  
  template <typename T1, typename A1, typename T2, typename A2>
  inline
  FeatureVectorLinear<T1,A1> operator*(FeatureVectorLinear<T1,A1>& x, const FeatureVector<T2,A2>& y)
  {
    typedef FeatureVectorLinear<T1,A1> self_type;
    
    self_type x_new(x);
    x_new *= y;
    return x_new;
  }

  template <typename T1, typename A1, typename T2, typename A2>
  inline
  FeatureVectorLinear<T1,A1> operator+(FeatureVectorLinear<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y)
  {
    typedef FeatureVectorLinear<T1,A1> self_type;
    
    self_type x_new(x);
    x_new += y;
    return x_new;
  }
  
  template <typename T1, typename A1, typename T2, typename A2>
  inline
  FeatureVectorLinear<T1,A1> operator-(FeatureVectorLinear<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y)
  {
    typedef FeatureVectorLinear<T1,A1> self_type;
    
    self_type x_new(x);
    x_new -= y;
    return x_new;    
  }
  
  template <typename T1, typename A1, typename T2, typename A2>
  inline
  FeatureVectorLinear<T1,A1> operator*(FeatureVectorLinear<T1,A1>& x, const FeatureVectorLinear<T2,A2>& y)
  {
    typedef FeatureVectorLinear<T1,A1> self_type;
    
    self_type x_new;
    self_type::multiply_equal_ordered(x_new.__map, x.begin(), x.end(), y.begin(), y.end());
    return x_new;
  }

};

namespace std
{
  template <typename T, typename A>
  inline
  void swap(cicada::FeatureVectorLinear<T,A>& x, cicada::FeatureVectorLinear<T, A>& y)
  {
    x.swap(y);
  }
};

#include <cicada/weight_vector.hpp>
#include <cicada/feature_vector.hpp>
#include <cicada/feature_vector_compact.hpp>

namespace cicada
{
  template <typename T, typename A>
  inline
  void FeatureVectorLinear<T,A>::assign(const FeatureVectorCompact& x)
  {
    __map.clear();
    __map.insert(x.begin(), x.end());
  }

  template <typename T, typename A>
  inline
  FeatureVectorLinear<T,A>& FeatureVectorLinear<T,A>::operator+=(const FeatureVectorCompact& x)
  {
    if (x.empty())
      return *this;
    else if (empty()) {
      assign(x);
      return *this;
    } else {
      plus_equal_ordered(__map, x.begin(), x.end());
      return *this;
    }
  }

  template <typename T, typename A>
  inline
  FeatureVectorLinear<T,A>& FeatureVectorLinear<T,A>::operator-=(const FeatureVectorCompact& x)
  {
    minus_equal_ordered(__map, x.begin(), x.end());
    return *this;
  }
  
  template <typename T, typename A>
  inline
  FeatureVectorLinear<T,A>& FeatureVectorLinear<T,A>::operator*=(const FeatureVectorCompact& x)
  {
    if (empty() || x.empty()) {
      clear();
      return *this;
    } else {
      map_type map_new;
      multiply_equal_ordered(map_new, __map.begin(), __map.end(), x.begin(), x.end());
      __map.swap(map_new);
      return *this;
    } 
  }
  
  template <typename T1, typename A1>
  inline
  FeatureVectorLinear<T1,A1> operator+(const FeatureVectorLinear<T1,A1>& x, const FeatureVectorCompact& y)
  {
    typedef FeatureVectorLinear<T1,A1> self_type;
    
    self_type x_new(x);
    x_new += y;
    return x_new;
  }
  
  template <typename T1, typename A1>
  inline
  FeatureVectorLinear<T1,A1> operator-(const FeatureVectorLinear<T1,A1>& x, const FeatureVectorCompact& y)
  {
    typedef FeatureVectorLinear<T1,A1> self_type;
    
    self_type x_new(x);
    x_new -= y;
    return x_new;
  }

  template <typename T1, typename A1>
  inline
  FeatureVectorLinear<T1,A1> operator*(const FeatureVectorLinear<T1,A1>& x, const FeatureVectorCompact& y)
  {
    typedef FeatureVectorLinear<T1,A1> self_type;
    
    self_type x_new;
    self_type::multiply_equal_ordered(x_new.__map, x.begin(), x.end(), y.begin(), y.end());
    return x_new;
  }

};

#endif
