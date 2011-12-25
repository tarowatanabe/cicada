// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
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

#include <google/dense_hash_map>

//
// feature vector, for use as "temporary" unordered vector for faster access, w/o sorting
//

namespace cicada
{
  
  // forward declaration...
  template <typename Tp, typename Alloc >
  class FeatureVector;
  
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
    FeatureVectorUnordered(const FeatureVectorUnordered<T,A>& x) : __map() { initialize(); __map.insert(x.begin(), x.end()); }
    template <typename Iterator>
    FeatureVectorUnordered(Iterator first, Iterator last) : __map() { initialize(); __map.insert(first, last); } 
     
    FeatureVectorUnordered& operator=(const self_type& x)
    {
      __map = x.__map;
      return *this;
    }
     
    template <typename T, typename A>
    FeatureVectorUnordered& operator=(const FeatureVectorUnordered<T,A>& x)
    {
      __map.clear();
      __map.insert(x.begin(), x.end());
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
    FeatureVectorUnordered& operator+=(const FeatureVectorUnordered<T,A>& x)
    {
      typedef FeatureVectorUnordered<T,A> another_type;
       
      typename another_type::const_iterator iter_end = x.end();
      for (typename another_type::const_iterator iter = x.begin(); iter != iter_end; ++ iter) {
	std::pair<iterator, bool> result = __map.insert(*iter);
	if (! result.second) {
	  result.first->second += iter->second;
	   
	  if (result.first->second == Tp())
	    __map.erase(result.first);
	}
      }
       
      return *this;
    }

    template <typename T, typename A>
    FeatureVectorUnordered& operator-=(const FeatureVectorUnordered<T,A>& x)
    {
      typedef FeatureVectorUnordered<T,A> another_type;
       
      typename another_type::const_iterator iter_end = x.end();
      for (typename another_type::const_iterator iter = x.begin(); iter != iter_end; ++ iter) {
	std::pair<iterator, bool> result = __map.insert(std::make_pair(iter->first, - iter->second));
	if (! result.second) {
	  result.first->second -= iter->second;
	   
	  if (result.first->second == Tp())
	    __map.erase(result.first);
	}
      }
       
      return *this;
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


#endif
