// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SEMIRING__LOG__HPP__
#define __CICADA__SEMIRING__LOG__HPP__ 1

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
    class Log
    {
    public:
      typedef Tp weight_type;
      typedef Tp value_type;
      typedef Log<Tp> self_type;
      
    public:
      struct proxy_type
      {
	friend struct Log;
	
	proxy_type(const weight_type& x, const bool& s) : __value(x), __sign(s) {}
	
	operator Log() const { return Log(*this); }
	
      private:
	weight_type __value;
	char        __sign;
      };

    public:
      static inline self_type log(const weight_type& x, const bool& s) { return proxy_type(x, s); }
      
      static inline self_type zero() { return proxy_type(impl::traits_infinity<value_type>::minus(), false); }
      static inline self_type one()  { return proxy_type(0, false); }
      static inline self_type max()  { return proxy_type(impl::traits_infinity<value_type>::plus(), false); }
      static inline self_type min()  { return proxy_type(impl::traits_infinity<value_type>::plus(), true); }
      
    public:
      Log() : __value(impl::traits_infinity<value_type>::minus()), __sign(false) {}
      Log(const weight_type& x) : __value(std::signbit(x) ? std::log(-x) : std::log(x)), __sign(std::signbit(x)) {}
      explicit Log(const proxy_type& x) : __value(x.__value), __sign(x.__sign) {}

      operator Tp() const { return  __sign ? - std::exp(__value) : std::exp(__value); }
      
    public:
      template <typename T>
      friend
      T log(const Log<T>& x);
      
      Log& operator+=(const Log& x)
      {
	using namespace boost::math::policies;
	typedef policy<domain_error<errno_on_error>,
	  pole_error<errno_on_error>,
	  overflow_error<errno_on_error>,
	  rounding_error<errno_on_error>,
	  evaluation_error<errno_on_error>
	  > policy_type;

	if (x.__value == impl::traits_infinity<value_type>::minus())
	  return *this;
	else if (__value == impl::traits_infinity<value_type>::minus()) {
	  *this = x;
	  return *this;
	}

	if (__sign == x.__sign) {
	  if (x.__value < __value)
	    __value = __value + boost::math::log1p(std::exp(x.__value - __value));
	  else
	    __value = x.__value + boost::math::log1p(std::exp(__value - x.__value));
	  
	} else {
	  if (x.__value == __value)
	    *this = zero();
	  else if (x.__value < __value) {
	    const Tp exp_value = std::exp(x.__value - __value);
	    if (exp_value == 1.0)
	      *this = zero();
	    else
	      __value = __value + boost::math::log1p(- exp_value);
	  } else {
	    const Tp exp_value = std::exp(__value - x.__value);
	    if (exp_value == 1.0)
	      *this = zero();
	    else {
	      __value = x.__value + boost::math::log1p(- exp_value);
	      __sign = ! __sign;
	    }
	  }
	}
	
	return *this;
      }
      
      Log& operator-=(const Log& x)
      {
	return *this += Log(proxy_type(x.__value, ! x.__sign));
      }
      
      Log& operator*=(const Log& x)
      {
	__sign = (__sign != x.__sign);
	__value += x.__value;
	return *this;
      }

      Log& operator/=(const Log& x)
      {
	__sign = (__sign != x.__sign);
	__value -= x.__value;
	return *this;
      }
      
      friend
      bool operator==(const self_type& x, const self_type& y) { return x.__value == y.__value && x.__sign == y.__sign; }
      friend
      bool operator!=(const self_type& x, const self_type& y) { return x.__value != y.__value || x.__sign != y.__sign; }
      friend
      bool operator<(const self_type& x, const self_type& y) { return (x.__sign > y.__sign) || (x.__sign && x.__value > y.__value) || (x.__value < y.__value); }
      friend
      bool operator>(const self_type& x, const self_type& y) { return y < x; }
      friend
      bool operator<=(const self_type& x, const self_type& y) { return ! (y < x); }
      friend
      bool operator>=(const self_type& x, const self_type& y) { return ! (x < y); }
      
      friend
      std::ostream& operator<<(std::ostream& os, const self_type& x)
      {
	os << (x.__sign ? '-' : '+') << x.__value;
	return os;
      }
      
      friend
      std::istream& operator>>(std::istream& is, self_type& x)
      {
	char __char;

	is >> __char >> x.__value;

	switch (__char) {
	case '+': x.__sign = false; break;
	case '-': x.__sign = true;  break;
	default: 
	  throw std::runtime_error("invlaid sign");
	}
	return is;
      }

      
    private:
      weight_type __value;
      char        __sign;
    };

    template <typename Tp>
    inline
    Tp log(const Log<Tp>& x)
    {
      if (x.__sign)
	throw std::runtime_error("no negative log");
      
      return x.__value;
    }
    
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
    Log<Tp> operator-(const Log<Tp>& x, const Log<Tp>& y)
    {
      Log<Tp> __value(x);
      __value -= y;
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
    inline
    Log<Tp> operator/(const Log<Tp>& x, const Log<Tp>& y)
    {
      Log<Tp> __value(x);
      __value /= y;
      return __value;
    }
    
    template <typename Tp>
    struct traits<Log<Tp> >
    {
      static inline Log<Tp> log(const Tp& x) { return Log<Tp>::log(x, false); }
      static inline Log<Tp> zero() { return Log<Tp>::zero();  }
      static inline Log<Tp> one()  { return Log<Tp>::one(); }
    };
  };
};


#endif
