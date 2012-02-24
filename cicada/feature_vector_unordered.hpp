// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE_VECTOR_UNORDERED__HPP__
#define __CICADA__FEATURE_VECTOR_UNORDERED__HPP__ 1

#include <map>
#include <memory>
#include <utility>
#include <algorithm>
#include <iterator>
#include <iostream>

#include <cicada/feature.hpp>

#include <utils/dense_hash_map.hpp>

//
// feature vector, for use as "temporary" unordered vector for faster access, w/o sorting
//

namespace cicada
{
  
  // forward declaration...
  template <typename Tp, typename Alloc >
  class FeatureVector;

  class FeatureVectorCompact;
  
  template <typename Tp, typename Alloc >
  class FeatureVectorLinear;
  
  template <typename Tp, typename Alloc >
  class WeightVector;
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class FeatureVectorUnordered
  {
  public:
    typedef cicada::Feature feature_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

  private:
    typedef google::dense_hash_map<feature_type, Tp, boost::hash<feature_type>, std::equal_to<feature_type> > map_type;
     
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
    typedef FeatureVectorUnordered<Tp, Alloc> self_type;
     
  public:
    FeatureVectorUnordered() : __map() { initialize(); }
    FeatureVectorUnordered(const self_type& x) : __map(x.__map) { }
    template <typename T, typename A>
    FeatureVectorUnordered(const FeatureVector<T,A>& x) : __map() { initialize(); assign(x); }
    FeatureVectorUnordered(const FeatureVectorCompact& x) : __map() { initialize(); assign(x); }
    template <typename T, typename A>
    FeatureVectorUnordered(const FeatureVectorLinear<T,A>& x) : __map() { initialize(); assign(x); }
    template <typename T, typename A>
    FeatureVectorUnordered(const FeatureVectorUnordered<T,A>& x) : __map() { initialize(); assign(x); }
    template <typename Iterator>
    FeatureVectorUnordered(Iterator first, Iterator last) : __map() { initialize(); __map.insert(first, last); } 
     
    FeatureVectorUnordered& operator=(const self_type& x)
    {
      __map = x.__map;
      return *this;
    }

    template <typename T, typename A>
    FeatureVectorUnordered& operator=(const FeatureVector<T,A>& x)
    {
      assign(x);
      return *this;
    }

    FeatureVectorUnordered& operator=(const FeatureVectorCompact& x)
    {
      assign(x);
      return *this;
    }
    
    template <typename T, typename A>
    FeatureVectorUnordered& operator=(const FeatureVectorLinear<T,A>& x)
    {
      assign(x);
      return *this;
    }

    template <typename T, typename A>
    FeatureVectorUnordered& operator=(const FeatureVectorUnordered<T,A>& x)
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
      
      if (x.sparse())
	__map.insert(x.sparse_begin(), x.sparse_end());
      else
	__map.insert(x.dense_begin(), x.dense_end());
    }

    void assign(const FeatureVectorCompact& x);
    
    template <typename T, typename A>
    void assign(const FeatureVectorLinear<T,A>& x)
    {
      __map.clear();
      __map.insert(x.begin(), x.end());
    }
     
    template <typename T, typename A>
    void assign(const FeatureVectorUnordered<T,A>& x)
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
     
    inline const_iterator find(const key_type& x) const { return __map.find(x); }
    inline       iterator find(const key_type& x)       { return __map.find(x); }

    void erase(const key_type& x) { __map.erase(x); }
     
    void swap(FeatureVectorUnordered& x) { __map.swap(x.__map); }
     
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
    FeatureVectorUnordered& operator+=(const T& x)
    {
      std::for_each(begin(), end(), __apply_unary<std::plus<Tp>, T>(x));
      return *this;
    }

    template <typename T>
    FeatureVectorUnordered& operator-=(const T& x)
    {
      std::for_each(begin(), end(), __apply_unary<std::minus<Tp>, T>(x));
      return *this;
    }
     
    template <typename T>
    FeatureVectorUnordered& operator*=(const T& x)
    {
      if (x == T())
	clear();
      else
	std::for_each(begin(), end(), __apply_unary<std::multiplies<Tp>, T>(x));
      return *this;
    }
     
    template <typename T>
    FeatureVectorUnordered& operator/=(const T& x)
    {
      std::for_each(begin(), end(), __apply_unary<std::divides<Tp>, T>(x));
      return *this;
    }
    

    template <typename T, typename A>
    FeatureVectorUnordered& operator+=(const FeatureVector<T,A>& x)
    {
      if (x.sparse())
	plus_equal(x.sparse_begin(), x.sparse_end());
      else
	plus_equal(x.dense_begin(), x.dense_end());
      return *this;
    }
    
    FeatureVectorUnordered& operator+=(const FeatureVectorCompact& x);

    template <typename T, typename A>
    FeatureVectorUnordered& operator+=(const FeatureVectorLinear<T,A>& x)
    {
      plus_equal(x.begin(), x.end());
      return *this;
    }
     
    template <typename T, typename A>
    FeatureVectorUnordered& operator+=(const FeatureVectorUnordered<T,A>& x)
    {
      plus_equal(x.begin(), x.end());
      return *this;
    }

    template <typename T, typename A>
    FeatureVectorUnordered& operator-=(const FeatureVector<T,A>& x)
    {
      if (x.sparse())
	minus_equal(x.sparse_begin(), x.sparse_end());
      else
	minus_equal(x.dense_begin(), x.dense_end());
      return *this;
    }

    FeatureVectorUnordered& operator-=(const FeatureVectorCompact& x);

    template <typename T, typename A>
    FeatureVectorUnordered& operator-=(const FeatureVectorLinear<T,A>& x)
    {
      minus_equal(x.begin(), x.end());
      return *this;
    }

    template <typename T, typename A>
    FeatureVectorUnordered& operator-=(const FeatureVectorUnordered<T,A>& x)
    {
      minus_equal(x.begin(), x.end());
      return *this;
    }

  public:
    friend bool operator==(const FeatureVectorUnordered& x, const FeatureVectorUnordered& y) { return x.__map == y.__map; }
    friend bool operator!=(const FeatureVectorUnordered& x, const FeatureVectorUnordered& y) { return x.__map != y.__map; }

  private:
    template <typename Iterator>
    void plus_equal(Iterator first, Iterator last)
    {
      for (/**/; first != last; ++ first) {
	std::pair<iterator, bool> result = __map.insert(*first);
	if (! result.second) {
	  result.first->second += first->second;
	  
	  if (result.first->second == Tp())
	    __map.erase(result.first);
	}
      }
    }
    
    template <typename Iterator>
    void minus_equal(Iterator first, Iterator last)
    {
      for (/**/; first != last; ++ first) {
	std::pair<iterator, bool> result = __map.insert(std::make_pair(first->first, - Tp(first->second)));
	if (! result.second) {
	  result.first->second -= first->second;
	  
	  if (result.first->second == Tp())
	    __map.erase(result.first);
	}
      }
    }

  private:
    void initialize()
    {
      __map.set_empty_key(feature_type(feature_type::id_type(-1)));
      __map.set_deleted_key(feature_type(feature_type::id_type(-2)));
    }
     
  private:
    map_type __map;
  };
  
  template <typename T1, typename A1, typename T2>
  inline
  FeatureVectorUnordered<T1,A1> operator+(const FeatureVectorUnordered<T1,A1>& x, const T2& y)
  {
    FeatureVectorUnordered<T1,A1> features(x);
    features += y;
    return features;
  }

  template <typename T2, typename T1, typename A1>
  inline
  FeatureVectorUnordered<T1,A1> operator+(const T2& x, const FeatureVectorUnordered<T1,A1>& y)
  {
    FeatureVectorUnordered<T1,A1> features(y);
    features += x;
    return features;
  }

  template <typename T1, typename A1, typename T2>
  inline
  FeatureVectorUnordered<T1,A1> operator*(const FeatureVectorUnordered<T1,A1>& x, const T2& y)
  {
    if (y == T2()) return FeatureVectorUnordered<T1,A1>();
    
    FeatureVectorUnordered<T1,A1> features(x);
    features *= y;
    return features;
  }

  template <typename T2, typename T1, typename A1>
  inline
  FeatureVectorUnordered<T1,A1> operator*(const T2& x, const FeatureVectorUnordered<T1,A1>& y)
  {
    if (x == T2()) return FeatureVectorUnordered<T1,A1>();
    
    FeatureVectorUnordered<T1,A1> features(y);
    features *= x;
    return features;
  }

  template <typename T1, typename A1, typename T2>
  inline
  FeatureVectorUnordered<T1,A1> operator-(const FeatureVectorUnordered<T1,A1>& x, const T2& y)
  {
    FeatureVectorUnordered<T1,A1> features(x);
    features -= y;
    return features;
  }
  
  template <typename T1, typename A1, typename T2>
  inline
  FeatureVectorUnordered<T1,A1> operator/(const FeatureVectorUnordered<T1,A1>& x, const T2& y)
  {
    FeatureVectorUnordered<T1,A1> features(x);
    features /= y;
    return features;
  }
};

namespace std
{
  template <typename T, typename A>
  inline
  void swap(cicada::FeatureVectorUnordered<T,A>& x, cicada::FeatureVectorUnordered<T, A>& y)
  {
    x.swap(y);
  }
};

#include <cicada/weight_vector.hpp>
#include <cicada/feature_vector.hpp>
#include <cicada/feature_vector_compact.hpp>
#include <cicada/feature_vector_linear.hpp>

namespace cicada
{
  template <typename T, typename A>
  inline
  void FeatureVectorUnordered<T,A>::assign(const FeatureVectorCompact& x)
  {
    __map.clear();
    __map.insert(x.begin(), x.end());
  }
  
  template <typename T, typename A>
  inline
  FeatureVectorUnordered<T,A>& FeatureVectorUnordered<T,A>::operator+=(const FeatureVectorCompact& x)
  {
    plus_equal(x.begin(), x.end());
    return *this;
  }

  template <typename T, typename A>
  inline
  FeatureVectorUnordered<T,A>& FeatureVectorUnordered<T,A>::operator-=(const FeatureVectorCompact& x)
  {
    minus_equal(x.begin(), x.end());
    return *this;
  }

};

#endif
