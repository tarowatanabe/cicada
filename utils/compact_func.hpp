// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__COMPACT_FUNC__HPP__
#define __UTILS__COMPACT_FUNC__HPP__

namespace utils
{
  template <typename Tp>
  struct unassigned
  {
    Tp operator()() { return Tp(); }
  };
  
  template <typename Tp>
  struct deleted
  {
    Tp operator()() { return Tp(); }
  };
  
};

#endif
