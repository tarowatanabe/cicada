// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__COMPACT_FUNC__HPP__
#define __UTILS__COMPACT_FUNC__HPP__

#include <utils/traits.hpp>

namespace utils
{
  template <typename Tp>
  struct unassigned
  {
    Tp operator()() const { return utils::traits<Tp>::unassigned(); }
  };
  
  template <typename Tp>
  struct deleted
  {
    Tp operator()() const { return utils::traits<Tp>::deleted(); }
  };
  
};

#endif
