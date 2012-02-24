// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SPARSE_HASH_MAP__HPP__
#define __UTILS__SPARSE_HASH_MAP__HPP__

#include <utils/config.hpp>

#include <boost/functional/hash.hpp>

#if defined(HAVE_SPARSEHASH)
  #include <sparsehash/sparse_hash_map>
#elif defined(HAVE_GOOGLE_SPARSEHASH)
  #include <google/sparse_hash_map>
#else
  #error "no sparse_hash_map?"
#endif

namespace utils
{
  template <typename _Key, typename _Tp,
	    typename _Hash=boost::hash<_Key>,
	    typename _Pred=std::equal_to<_Key>,
	    typename _Alloc=std::allocator<std::pair<const _Key, _Tp> > >
  struct sparse_hash_map
  {
    typedef google::sparse_hash_map<_Key,_Tp,_Hash,_Pred> type;
  };
};

#endif
