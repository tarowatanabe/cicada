// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SPARSE_HASH_MAP__HPP__
#define __UTILS__SPARSE_HASH_MAP__HPP__

#include <utils/config.hpp>

#if defined(HAVE_SPARSEHASH)
  #include <sparsehash/sparse_hash_map>
#elif defined(HAVE_GOOGLE_SPARSEHASH)
  #include <google/sparse_hash_map>
#else
  #error "no sparse_hash_map?"
#endif

#endif
