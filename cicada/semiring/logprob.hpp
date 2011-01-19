// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SEMIRING__LOGPROB__HPP__
#define __CICADA__SEMIRING__LOGPROB__HPP__ 1

#include <cmath>
#include <cfloat>
#include <climits>

#include <limits>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <boost/math/special_functions/log1p.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include <cicada/semiring/traits.hpp>

namespace cicada
{
  namespace semiring
  {

    template <typename Tp>
    class Logprob
    {
    public:
      typedef Tp weight_type;
      typedef Tp value_type;
      typedef Logprob<Tp> self_type;

    public:
      struct proxy_type
      {
	friend struct Logprob;
	
	proxy_type(const weight_type& x) : __value(x) {}
	
	operator Logprob() const { return Logprob(*this); }
	
      private:
	weight_type __value;
      };
      
    public:
      static inline self_type exp(const Tp& x) { return proxy_type(x); }
      
      static inline self_type zero() { return proxy_type(impl::traits_infinity<value_type>::minus()); }
      static inline self_type one()  { return proxy_type(0); }
      static inline self_type max()  { return proxy_type(impl::traits_infinity<value_type>::plus()); }
      static inline self_type min()  { return proxy_type(impl::traits_infinity<value_type>::minus()); }
      
    public:
      // any better wayt to make an assignment...?
      Logprob() : __value(impl::traits_infinity<value_type>::minus()) {}
      Logprob(const weight_type& x) : __value(std::log(x)) {}
      explicit Logprob(const proxy_type& x) : __value(x.__value) {}


      operator Tp() const { return std::exp(__value); }
      
    public:
      template <typename T>
      friend
      T log(const Logprob<T>& x);

      Logprob& operator+=(const Logprob& x)
      {
	if (*this == zero()) {
	  __value = x.__value;
	  return *this;
	} else if (x == zero())
	  return *this;
	
	if (__value <= x.__value) 
	  __value = x.__value + boost::math::log1p(std::exp(__value - x.__value));
	else
	  __value = __value + boost::math::log1p(std::exp(x.__value - __value));
	
	return *this;
      }

      Logprob& operator-=(const Logprob& x)
      {
	if (x > *this)
	  throw std::runtime_error("invalid minus");
	
	if (*this == x)
	  *this = zero();
	else {
	  const Tp exp_value = std::exp(x.__value - __value);
	  if (exp_value == 1.0)
	    *this = zero();
	  else
	    __value = __value + boost::math::log1p(- exp_value);
	}
	
	return *this;
      }
      
      Logprob& operator*=(const Logprob& x)
      {
	__value += x.__value;
	return *this;
      }

      Logprob& operator/=(const Logprob& x)
      {
	__value -= x.__value;
	return *this;
      }
      
      friend
      bool operator==(const self_type& x, const self_type& y) { return x.__value == y.__value; }
      friend
      bool operator!=(const self_type& x, const self_type& y) { return x.__value != y.__value; }
      friend
      bool operator>(const self_type& x, const self_type& y) { return x.__value > y.__value; }
      friend
      bool operator<(const self_type& x, const self_type& y) { return x.__value < y.__value; }
      friend
      bool operator>=(const self_type& x, const self_type& y) { return x.__value >= y.__value; }
      friend
      bool operator<=(const self_type& x, const self_type& y) { return x.__value <= y.__value; }

      friend
      std::ostream& operator<<(std::ostream& os, const self_type& x)
      {
	os << x.__value;
	return os;
      }
      
      friend
      std::istream& operator>>(std::istream& is, self_type& x)
      {
	is >> x.__value;
	return is;
      }
      
    private:
      weight_type __value;
    };

    template <typename Tp>
    inline
    Tp log(const Logprob<Tp>& x)
    {
      return x.__value;
    }
    
    template <typename Tp>
    inline
    Logprob<Tp> operator+(const Logprob<Tp>& x, const Logprob<Tp>& y)
    {
      Logprob<Tp> __value(x);
      __value += y;
      return __value;
    }

    template <typename Tp>
    inline
    Logprob<Tp> operator-(const Logprob<Tp>& x, const Logprob<Tp>& y)
    {
      Logprob<Tp> __value(x);
      __value -= y;
      return __value;
    }

    template <typename Tp>
    inline
    Logprob<Tp> operator*(const Logprob<Tp>& x, const Logprob<Tp>& y)
    {
      Logprob<Tp> __value(x);
      __value *= y;
      return __value;
    }

    template <typename Tp>
    inline
    Logprob<Tp> operator/(const Logprob<Tp>& x, const Logprob<Tp>& y)
    {
      Logprob<Tp> __value(x);
      __value /= y;
      return __value;
    }
    
    template <typename Tp>
    struct traits<Logprob<Tp> >
    {
      static inline Logprob<Tp> exp(const Tp& x) { return Logprob<Tp>::exp(x); }
      static inline Logprob<Tp> zero() { return Logprob<Tp>::zero();  }
      static inline Logprob<Tp> one()  { return Logprob<Tp>::one(); }
    };

  };
};


#endif
