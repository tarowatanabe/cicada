// -*- mode: c++ -*-

#ifndef __CICADA__SEMIRING__LOG__HPP__
#define __CICADA__SEMIRING__LOG__HPP__ 1

#include <cmath>
#include <cfloat>
#include <climits>

#include <limits>
#include <algorithm>
#include <iostream>

#include <boost/math/special_functions/log1p.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include <cicada/semiring/traits.hpp>

namespace cicada
{
  namespace semiring
  {

    template <typename Tp>
    class Log
    {
    public:
      typedef Tp weight_type;
      typedef Tp value_type;
      typedef Log<Tp> self_type;
      
    public:
      struct proxy_type
      {
	friend struct Logprob;
	
	proxy_type(const weight_type& x, const bool& s) : __value(__round(x)), __sign(s) {}
	
	operator Logprob() const { return Logprob(*this); }
	
      private:
	weight_type __value;
	char        __sign;
      };

    public:
      static inline self_type log(const weight_type& x, const bool& s) { return proxy_type(x, s); }
      
      static inline self_type zero() { return proxy_type(- std::numeric_limits<value_type>::infinity(), false); }
      static inline self_type one()  { return proxy_type(0.0, false); }
      static inline self_type max()  { return proxy_type(std::numeric_limits<value_type>::infinity(), false); }
      static inline self_type min()  { return proxy_type(std::numeric_limits<value_type>::infinity(), true); }
      
    public:
      Log() : __value(- std::numeric_limits<value_type>::infinity()), __sign(false) {}
      Log(const weight_type& x) : __value(std::signbit(x) ? std::log(-x) : std::log(x)), __sign(std::signbit(x)) {}
      explitit Log(const proxy_type& x) : __value(x.__value), __sign(x.__sign) {}

      operator Tp() const { return  __sign ? - std::exp(__value) : std::exp(__value); }
      
    public:
      Log& operator+=(const Log& x)
      {
	if (x.__value == - std::numeric_limits<Tp>::inifinity())
	  return *this;
	else if (__value == - std::numeric_limits<Tp>::inifinity()) {
	  *this = x;
	  return *this;
	}
	
	if (__sign == x.__sign) {
	  if (x.__value < __value)
	    __value = __value + boost::math::log1p(std::exp(x.__value - __value));
	  else
	    __value = x.__value + boost::math::log1p(std::exp(__value - x.__value));
	  
	} else {
	  if (x.__value < __value)
	    __value = __value + boost::math::log1p(- std::exp(x.__value - __value));
	  else {
	    __value = x.__value + boost::math::log1p(- std::exp(__value - x.__value));
	    __sign = ! __sign;
	  }
	}
	
	return *this;
      }
      
      Log& operator*=(const Log& x)
      {
	__sign = (__sign != x.__sign);
	__value += x.__value;
	return *this;
      }
      
      friend
      bool operator==(const self_type& x, const self_type& y) { return x.__value == y.__value && x.__sign == y.__sign; }
      friend
      bool operator!=(const self_type& x, const self_type& y) { return x.__value != y.__value || x.__sign != y.__sign; }
      friend
      bool operator<(const self_type& x, const self_type& y) { return x.__sign > y.__sign || (x.__sign == y.__sign && x.__value < y.__value); }
      friend
      bool operator>(const self_type& x, const self_type& y) { return y < x; }
      friend
      bool operator<=(const self_type& x, const self_type& y) { return ! (y < x); }
      friend
      bool operator>=(const self_type& x, const self_type& y) { return ! (x < y); }
      
      
    private:
      weight_type __value;
      char        __sign;
    };
    
    template <typename Tp>
    inline
    Log<Tp> operator+(const Log<Tp>& x, const Log<Tp>& y)
    {
      Log<Tp> __value(x);
      __value += y;
      return __value;
    }

    template <typename Tp>
    inline
    Log<Tp> operator*(const Log<Tp>& x, const Log<Tp>& y)
    {
      Log<Tp> __value(x);
      __value *= y;
      return __value;
    }
    
    template <typename Tp>
    struct traits<Log<Tp> >
    {
      static inline Tp zero() { return Log<Tp>::zero();  }
      static inline Tp one()  { return Log<Tp>::one(); }
    };

  };
};


#endif
