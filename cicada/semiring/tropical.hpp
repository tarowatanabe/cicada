// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SEMIRING__TROPICAL__HPP__
#define __CICADA__SEMIRING__TROPICAL__HPP__ 1

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
    class Tropical
    {
    public:
      typedef Tp weight_type;
      typedef Tp value_type;
      typedef Tropical<Tp> self_type;

    public:
      struct proxy_type
      {
	friend struct Tropical;
	
	proxy_type(const weight_type& x) : __value(x) {}
	
	operator Tropical() const { return Tropical(*this); }
	
      private:
	weight_type __value;
      };
      
    public:
      static inline self_type log(const Tp& x) { return proxy_type(x); }
      
      static inline self_type zero() { return proxy_type(impl::traits_infinity<value_type>::minus()); }
      static inline self_type one()  { return proxy_type(0); }
      static inline self_type max()  { return proxy_type(impl::traits_infinity<value_type>::plus()); }
      static inline self_type min()  { return proxy_type(impl::traits_infinity<value_type>::minus()); }
      
    public:
      // any better wayt to make an assignment...?
      Tropical() : __value(impl::traits_infinity<value_type>::minus()) {}
      Tropical(const weight_type& x) : __value(std::log(x)) {}
      explicit Tropical(const proxy_type& x) : __value(x.__value) {}


      operator Tp() const { return std::exp(__value); }
      
    public:
      template <typename T>
      friend
      T log(const Tropical<T>& x);

      Tropical& operator+=(const Tropical& x)
      {
	__value = std::max(__value, x.__value);
	
	return *this;
      }
      
      Tropical& operator*=(const Tropical& x)
      {
	__value += x.__value;
	return *this;
      }

      Tropical& operator/=(const Tropical& x)
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
    Tp log(const Tropical<Tp>& x)
    {
      return x.__value;
    }
    
    template <typename Tp>
    inline
    Tropical<Tp> operator+(const Tropical<Tp>& x, const Tropical<Tp>& y)
    {
      Tropical<Tp> __value(x);
      __value += y;
      return __value;
    }

    template <typename Tp>
    inline
    Tropical<Tp> operator*(const Tropical<Tp>& x, const Tropical<Tp>& y)
    {
      Tropical<Tp> __value(x);
      __value *= y;
      return __value;
    }

    template <typename Tp>
    inline
    Tropical<Tp> operator/(const Tropical<Tp>& x, const Tropical<Tp>& y)
    {
      Tropical<Tp> __value(x);
      __value /= y;
      return __value;
    }
    
    template <typename Tp>
    struct traits<Tropical<Tp> >
    {
      static inline Tropical<Tp> log(const Tp& x) { return Tropical<Tp>::log(x); }
      static inline Tropical<Tp> zero() { return Tropical<Tp>::zero();  }
      static inline Tropical<Tp> one()  { return Tropical<Tp>::one(); }
    };

  };
};


#endif
