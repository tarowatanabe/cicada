// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SEMIRING__TRAITS__HPP__
#define __CICADA__SEMIRING__TRAITS__HPP__ 1

#include <cmath>

#include <boost/numeric/conversion/bounds.hpp>

namespace cicada
{
  namespace semiring
  {
    namespace impl
    {
      template <typename Tp, bool has_infinity>
      struct __traits_infinity {};
      
      template <typename Tp>
      struct __traits_infinity<Tp, true>
      {
	static const Tp plus()  { return std::numeric_limits<Tp>::infinity(); }
	static const Tp minus() { return - std::numeric_limits<Tp>::infinity(); }
      };
      
      template <typename Tp>
      struct __traits_infinity<Tp, false>
      {
	static const Tp plus()  { return std::numeric_limits<Tp>::max(); }
	static const Tp minus() { return std::numeric_limits<Tp>::min(); }
      };

      template <typename Tp>
      struct traits_infinity : public __traits_infinity<Tp, std::numeric_limits<Tp>::has_infinity>
      {
	
	
      };
    };

    template <typename Tp>
    struct traits
    {
      template <typename T>
      static inline Tp log(const T& x)  { return std::log(x); }
      static inline Tp zero() { return Tp();  }
      static inline Tp one()  { return Tp(1); }
    };

  };
};


#endif
