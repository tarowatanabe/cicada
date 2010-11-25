// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SEMIRING__EXPECTATION__HPP__
#define __CICADA__SEMIRING__EXPECTATION__HPP__ 1

#include <cmath>
#include <cfloat>
#include <climits>

#include <limits>
#include <algorithm>
#include <iostream>

#include <cicada/semiring/traits.hpp>


namespace cicada
{
  namespace semiring
  {

    template <typename Tp, typename Tr>
    struct Expectation
    {
    public:
      typedef Tp p_type;
      typedef Tr r_type;
      
      typedef Expectation<Tp, Tr> self_type;
      
    public:
      static inline self_type zero() { return self_type(); }
      static inline self_type one()  { return self_type(traits<p_type>::one()); }
      static inline self_type max()  { return self_type(traits<p_type>::max()); }
      static inline self_type min()  { return self_type(traits<p_type>::min()); }
      
      Expectation() : p(), r() {}
      explicit Expectation(const p_type& _p) : p(_p), r() {}
      Expectation(const p_type& _p, const r_type& _r) : p(_p), r(_r) {}
      
      Expectation& operator+=(const Expectation& x)
      {
	p += x.p;
	r += x.r;
	return *this;
      }
      
      Expectation& operator*=(const Expectation& x)
      {
	r = x.r * p + r * x.p;
	p *= p;
	return *this;
      }
      
    public:
      p_type p;
      r_type r;
    };
    
    template <typename Tp, typename Tr>
    inline
    Expectation<Tp, Tr> operator+(const Expectation<Tp, Tr>& x, const Expectation<Tp, Tr>& y)
    {
      Expectation<Tp, Tr> __value(x);
      __value += y;
      return __value;
    }
    
    template <typename Tp, typename Tr>
    inline
    Expectation<Tp, Tr> operator*(const Expectation<Tp, Tr>& x, const Expectation<Tp, Tr>& y)
    {
      Expectation<Tp, Tr> __value(x);
      __value *= y;
      return __value;
    }
    
    template <typename Tp, typename Tr>
    struct traits<Expectation<Tp, Tr> >
    {
      static inline Tp zero() { return Expectation<Tp, Tr>::zero();  }
      static inline Tp one()  { return Expectation<Tp, Tr>::one(); }
    };

  };
};

#endif
