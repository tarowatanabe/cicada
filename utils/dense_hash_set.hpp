// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__DENSE_HASH_SET__HPP__
#define __UTILS__DENSE_HASH_SET__HPP__

#include <utils/config.hpp>

#if defined(HAVE_SPARSEHASH)
  #include <sparsehash/dense_hash_set>
#elif defined(HAVE_GOOGLE_SPARSEHASH)
  #include <google/dense_hash_set>
#else
  #error "no dense_hash_set?"
#endif

#endif
