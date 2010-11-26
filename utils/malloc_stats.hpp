// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MALLOC_STATS__HPP__
#define __UTILS__MALLOC_STATS__HPP__ 1

#include <stdint.h>

#include <memory>

namespace utils
{
  struct malloc_stats
  {
    static size_t allocated(); // allocated memory
    static size_t used();      // actually used memory
  };
};

#endif
