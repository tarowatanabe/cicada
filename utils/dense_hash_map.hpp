// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__DENSE_HASH_MAP__HPP__
#define __UTILS__DENSE_HASH_MAP__HPP__

#include <utils/config.hpp>

#if defined(HAVE_SPARSEHASH)
  #include <sparsehash/dense_hash_map>
#elif defined(HAVE_GOOGLE_SPARSEHASH)
  #include <google/dense_hash_map>
#else
  #error "no dense_hash_map?"
#endif

#endif
