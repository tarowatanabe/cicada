// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__UNORDERED_MAP__HPP__
#define __UTILS__UNORDERED_MAP__HPP__ 1

#include <utils/config.hpp>

#include <boost/functional/hash.hpp>

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

namespace utils
{
  template <typename _Key, typename _Tp,
	    typename _Hash=boost::hash<_Key>,
	    typename _Pred=std::equal_to<_Key>,
	    typename _Alloc=std::allocator<std::pair<const _Key, _Tp> > >
  struct unordered_map
  {
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<_Key,_Tp,_Hash,_Pred,_Alloc> type;
#else
    typedef sgi::hash_map<_Key,_Tp,_Hash,_Pred,_Alloc> type;
#endif
  };

  template <typename _Key, typename _Tp,
	    typename _Hash=boost::hash<_Key>,
	    typename _Pred=std::equal_to<_Key>,
	    typename _Alloc=std::allocator<std::pair<const _Key, _Tp> > >
  struct unordered_multimap
  {
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_multimap<_Key,_Tp,_Hash,_Pred,_Alloc> type;
#else
    typedef sgi::hash_multimap<_Key,_Tp,_Hash,_Pred,_Alloc> type;
#endif
  };
};

#endif
