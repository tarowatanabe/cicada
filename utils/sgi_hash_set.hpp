// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SGI_HASH_SET__HPP__
#define __UTILS__SGI_HASH_SET__HPP__

#include <utils/config.hpp>

#if defined(HAVE_TR1_UNORDERED_SET)
  #if defined(__GNUC__) && ( (__GNUC__ > 4) || ((__GNUC__ >= 4) && (__GNUC_MINOR__ >= 1)) )
    #include <tr1/unordered_set>
  #else
    #undef HAVE_TR1_UNORDERED_SET
  #endif
#endif

#ifndef HAVE_TR1_UNORDERED_SET
  #if defined(HAVE_EXT_HASH_SET)
    #include <ext/hash_set>
    #if __GNUC__ == 3 && __GNUC_MINOR__ == 0
      namespace sgi = std;
    #else
      namespace sgi = ::__gnu_cxx;
    #endif
  #elif defined(HAVE_HASH_SET)
      #include <hash_set>
      namespace sgi
      {
	using ::hash_set;
	using ::hash_multiset;
      };
  #endif
#endif

#endif
