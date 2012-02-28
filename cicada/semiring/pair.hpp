// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SEMIRING__PAIR__HPP__
#define __CICADA__SEMIRING__PAIR__HPP__ 1

#include <limits>
#include <utility>

#include <cicada/semiring/traits.hpp>

namespace cicada
{
  namespace semiring
  {
    template <typename _First, typename _Second >
    class Pair
    {
    public:
      typedef _First  first_type;
      typedef _Second second_type;
      
    private:
      typedef Pair<_First, _Second> self_type;
      
    public:
      static inline self_type zero() { return self_type(); }
      static inline self_type one()  { return self_type(traits<first_type>::one(), traits<second_type>::one()); }
      
    public:
      Pair() : first(), second() {}
      Pair(const first_type& __first, const second_type& __second) : first(__first), second(__second) {}
      
    public:
      template <typename F, typename S>
      self_type& operator+=(const Pair<F,S>& x)
      {
	first  += x.first;
	second += x.second;
	
	return *this;
      }

      template <typename F, typename S>
      self_type& operator*=(const Pair<F,S>& x)
      {
	first  *= x.first;
	second *= x.second;
	
	return *this;
      }

      template <typename T>
      self_type& operator*=(const T& x)
      {
	first  *= x;
	second *= x;
	
	return *this;
      }
      
    public:
      first_type  first;
      second_type second;
    };
    
    
    template <typename _First, typename _Second, typename F, typename S>
    inline
    Pair<_First, _Second> operator+(const Pair<_First, _Second>& x, const Pair<F,S>& y)
    {
      Pair<_First, _Second> __value(x);
      __value += y;
      return __value;
    }
    
    template <typename _First, typename _Second, typename F, typename S>
    inline
    Pair<_First, _Second> operator*(const Pair<_First, _Second>& x, const Pair<F, S>& y)
    {
      Pair<_First, _Second> __value(x);
      __value *= y;
      return __value;
    }

    template <typename _First, typename _Second, typename T>
    inline
    Pair<_First, _Second> operator*(const Pair<_First, _Second>& x, const T& y)
    {
      Pair<_First, _Second> __value(x);
      __value *= y;
      return __value;
    }


    template <typename _First, typename _Second>
    struct traits<Pair<_First, _Second> >
    {
      static inline Pair<_First,_Second> zero() { return Pair<_First,_Second>::zero();  }
      static inline Pair<_First,_Second> one() { return Pair<_First,_Second>::one();  }
    };

  };
};

#endif
