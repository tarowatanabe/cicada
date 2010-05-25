// -*- mode: c++ -*-

#ifndef __CICADA__WEIGHT_VECTOR__HPP__
#define __CICADA__WEIGHT_VECTOR__HPP__ 1

#include <vector>
#include <utility>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <string>

#include <cicada/feature.hpp>

#include <utils/space_separator.hpp>
#include <utils/bithack.hpp>

#include <boost/tokenizer.hpp>
#include <boost/fusion/tuple.hpp>

namespace cicada
{
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class WeightVector
  {
  public:
    typedef cicada::Feature feature_type;
    typedef cicada::Feature key_type;

  private:
    typedef std::vector<Tp, Alloc> weight_vector_type;
    typedef WeightVector<Tp, Alloc> self_type;
    
  public:
    typedef typename weight_vector_type::size_type       size_type;
    typedef typename weight_vector_type::difference_type difference_type;
    
    typedef typename weight_vector_type::value_type      value_type;
    
    typedef typename weight_vector_type::const_iterator  const_iterator;
    typedef typename weight_vector_type::iterator        iterator;
    typedef typename weight_vector_type::const_reference const_reference;
    typedef typename weight_vector_type::reference       reference;

  public:
    WeightVector() {}
    WeightVector(const WeightVector& x) : __values(x.__values) {}
    
  public:
    size_type size() const { return __values.size(); }
    bool empty() const { return __values.empty(); }
    
    void reserve(size_type x) { __values.reserve(x); }
    void clear() { __values.clear(); }
    
    void swap(WeightVector& x) { __values.swap(x.__values); }
    
    Tp operator[](const key_type& x) const
    {
      return (x.id() < __values.size() ? __values[x.id()] : Tp());
    }
    
    Tp& operator[](const key_type& x) 
    {
      if (x.id() >= __values.size())
	__values.resize(x.id() + 1, Tp());
      return __values[x.id()];
    }
        
    const_iterator begin() const { return __values.begin(); }
    iterator begin() { return __values.begin(); }
    const_iterator end() const { return __values.end(); }
    iterator end() { return __values.end(); }
    
  public:
    // operators...
    template <typename T>
    self_type& operator+=(const T& x)
    {
      std::transform(begin(), end(), begin(), std::bind2nd(std::plus<Tp>(), x));
    }
    
    template <typename T>
    self_type& operator-=(const T& x)
    {
      std::transform(begin(), end(), begin(), std::bind2nd(std::minus<Tp>(), x));
    }
    
    template <typename T>
    self_type& operator*=(const T& x)
    {
      std::transform(begin(), end(), begin(), std::bind2nd(std::multiplies<Tp>(), x));
    }
    
    template <typename T>
    self_type& operator/=(const T& x)
    {
      std::transform(begin(), end(), begin(), std::bind2nd(std::divides<Tp>(), x));
    }
    
    template <typename T, typename A>
    self_type& operator+=(const WeightVector<T, A>& x)
    {
      if (size() < x.size()) {
	__values.reserve(x.size());
	__values.resize(x.size());
      }
      
      std::transform(begin(), begin() + utils::bithack::min(size(), x.size()), x.begin(), begin(), std::plus<Tp>());
    }
    
    template <typename T, typename A>
    self_type& operator-=(const WeightVector<T, A>& x)
    {
      if (size() < x.size()) {
	__values.reserve(x.size());
	__values.resize(x.size());
      }
      
      std::transform(begin(), begin() + utils::bithack::min(size(), x.size()), x.begin(), begin(), std::minus<Tp>());
    }
    
    template <typename T, typename A>
    self_type& operator*=(const WeightVector<T, A>& x)
    {
      if (size() < x.size()) {
	__values.reserve(x.size());
	__values.resize(x.size());
      }
      
      std::transform(begin(), begin() + utils::bithack::min(size(), x.size()), x.begin(), begin(), std::multiplies<Tp>());
    }
    
    template <typename T, typename A>
    self_type& operator/=(const WeightVector<T, A>& x)
    {
      if (size() < x.size()) {
	__values.reserve(x.size());
	__values.resize(x.size());
      }
      
      std::transform(begin(), begin() + utils::bithack::min(size(), x.size()), x.begin(), begin(), std::divides<Tp>());
    }
    
  public:
    
    template <typename T, typename A>
    friend
    std::ostream& operator<<(std::ostream& os, const WeightVector<T, A>& x)
    {
      typename WeightVector<T,A>::const_iterator iter_begin = x.begin();
      typename WeightVector<T,A>::const_iterator iter_end   = x.end();

      static const cicada::Feature __empty;
      
      for (typename WeightVector<T,A>::const_iterator iter = iter_begin; iter != iter_end; ++ iter)
	if (*iter != 0.0) {
	  cicada::Feature feature(iter - iter_begin);
	  if (feature != __empty)
	    os << feature << ' ' << *iter << '\n';
	}
      
      return os;
    }
    
    template <typename T, typename A>
    friend
    std::istream& operator>>(std::istream& is, WeightVector<T,A>& x)
    {
      typedef boost::tokenizer<utils::space_separator> tokenizer_type;
      
      std::string feature;
      T value;
      while ((is >> feature) && (is >> value))
	x[cicada::Feature(feature).id()] = value;
    }
    
  private:
    weight_vector_type __values;
  };
    

};

namespace std
{
  
  template <typename T, typename A>
  inline
  void swap(cicada::WeightVector<T,A>& x, cicada::WeightVector<T,A>& y)
  {
    x.swap(y);
  }
  
};

#endif
