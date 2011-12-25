// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
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
  class FeatureVectorUnordered;
  
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
    template <typename T, typename A>
    FeatureVectorLinear(const FeatureVectorUnordered<T,A>& x) : __map() { assign(x); }
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

    template <typename T, typename A>
    FeatureVectorLinear& operator=(const FeatureVectorUnordered<T,A>& x)
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
	__map.insert(x.sparse_begin(), x.saprse_end());
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
    
     
  private:
    map_type __map;
  };
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
#include <cicada/feature_vector_unordered.hpp>

namespace cicada
{
  template <typename T, typename A>
  inline
  void FeatureVectorLinear<T,A>::assign(const FeatureVectorCompact& x)
  {
    __map.clear();
    __map.insert(x.begin(), x.end());
  }
};

#endif
