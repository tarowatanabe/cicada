// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__DENSE_HASH_MAP__HPP__
#define __UTILS__DENSE_HASH_MAP__HPP__

#include <utils/config.hpp>

#include <boost/functional/hash.hpp>

#if defined(HAVE_SPARSEHASH)
  #include <sparsehash/dense_hash_map>
#elif defined(HAVE_GOOGLE_SPARSEHASH)
  #include <google/dense_hash_map>
#else
  #error "no dense_hash_map?"
#endif

namespace utils
{
  template <typename _Key, typename _Tp,
	    typename _Hash=boost::hash<_Key>,
	    typename _Pred=std::equal_to<_Key>,
	    typename _Alloc=std::allocator<std::pair<const _Key, _Tp> > >
  struct dense_hash_map
  {
    typedef google::dense_hash_map<_Key,_Tp,_Hash,_Pred> type;
  };
};

#endif
