// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// traits, specialied for floating points...

#ifndef __UTILS__FLOAT_TRAITS__HPP__
#define __UTILS__FLOAT_TRAITS__HPP__ 1

#include <limits>

#include <boost/numeric/conversion/bounds.hpp>

#include <utils/logprob.hpp>

namespace utils
{
  
  template <typename Tp>
  struct float_traits
  {
    static inline Tp zero() { return Tp(); } 
    static inline Tp max() { return boost::numeric::bounds<Tp>::highest(); }
    static inline Tp min() { return boost::numeric::bounds<Tp>::lowest(); }
  };
  
};

#endif

