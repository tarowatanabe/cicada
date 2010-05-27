// -*- mode: c++ -*-

#ifndef __UTILS__SGI_HASH_MAP__HPP__
#define __UTILS__SGI_HASH_MAP__HPP__

#include <utils/config.hpp>

#if defined(HAVE_TR1_UNORDERED_MAP)
  #if defined(__GNUC__) && ( (__GNUC__ > 4) || ((__GNUC__ >= 4) && (__GNUC_MINOR__ >= 1)) )
    #include <tr1/unordered_map>
  #else
    #undef HAVE_TR1_UNORDERED_MAP
  #endif
#endif

#if ! defined(HAVE_TR1_UNORDERED_MAP)
  #if defined(HAVE_EXT_HASH_MAP)
    #include <ext/hash_map>
    #if __GNUC__ == 3 && __GNUC_MINOR__ == 0
      namespace sgi = std;
    #else
      namespace sgi = ::__gnu_cxx;
    #endif
  #elif defined(HAVE_HASH_MAP)
      #include <hash_map>
      namespace sgi
      {
	using ::hash_map;
	using ::hash_multimap;
      };
  #endif
#endif

#endif
