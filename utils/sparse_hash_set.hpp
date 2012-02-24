// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SPARSE_HASH_SET__HPP__
#define __UTILS__SPARSE_HASH_SET__HPP__

#include <utils/config.hpp>

#include <boost/functional/hash.hpp>

#if defined(HAVE_SPARSEHASH)
  #include <sparsehash/sparse_hash_set>
#elif defined(HAVE_GOOGLE_SPARSEHASH)
  #include <google/sparse_hash_set>
#else
  #error "no sparse_hash_set?"
#endif

namespace utils
{
  template <typename _Value,
	    typename _Hash=boost::hash<_Value>,
	    typename _Pred=std::equal_to<_Value>,
	    typename _Alloc=std::allocator<_Value > >
  struct sparse_hash_set
  {
    typedef google::sparse_hash_set<_Value,_Hash,_Pred> type;
  };
};

#endif
