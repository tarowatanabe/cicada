// -*- mode: c++ -*-

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
