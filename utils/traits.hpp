// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__TRAITS__HPP__
#define __UTILS__TRAITS__HPP__

namespace utils
{
  template <typename Tp>
  struct traits
  {
    static inline Tp unassigned() { return Tp(); }
    static inline Tp deleted() { return Tp(); }
  };
};

#endif
