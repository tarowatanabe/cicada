// -*- mode: c++ -*-
//
//  Copyright(C) 2012-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__UNORDERED_SET__HPP__
#define __UTILS__UNORDERED_SET__HPP__ 1

#include <utils/config.hpp>

#include <boost/functional/hash.hpp>

#if defined(HAVE_UNORDERED_SET)
  #include <unordered_set>
#elif defined(HAVE_TR1_UNORDERED_SET)
  #if defined(__GNUC__) && ( (__GNUC__ > 4) || ((__GNUC__ >= 4) && (__GNUC_MINOR__ >= 1)) )
    #include <tr1/unordered_set>
  #else
    #undef HAVE_TR1_UNORDERED_SET
  #endif
#endif

#if ! defined(HAVE_TR1_UNORDERED_SET) && ! defined(HAVE_UNORDERED_SET)
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

namespace utils
{
  template <typename _Value,
	    typename _Hash=boost::hash<_Value>,
	    typename _Pred=std::equal_to<_Value>,
	    typename _Alloc=std::allocator<_Value> >
  struct unordered_set
  {
#if defined(HAVE_UNORDERED_SET)
    typedef std::unordered_set<_Value,_Hash,_Pred,_Alloc> type;
#elif defined(HAVE_TR1_UNORDERED_SET)
    typedef std::tr1::unordered_set<_Value,_Hash,_Pred,_Alloc> type;
#else
    typedef sgi::hash_set<_Value,_Hash,_Pred,_Alloc> type;
#endif
  };

  template <typename _Value,
	    typename _Hash=boost::hash<_Value>,
	    typename _Pred=std::equal_to<_Value>,
	    typename _Alloc=std::allocator<_Value> >
  struct unordered_multiset
  {
#if defined(HAVE_UNORDERED_SET)
    typedef std::unordered_multiset<_Value,_Hash,_Pred,_Alloc> type;
#elif defined(HAVE_TR1_UNORDERED_SET)
    typedef std::tr1::unordered_multiset<_Value,_Hash,_Pred,_Alloc> type;
#else
    typedef sgi::hash_multiset<_Value,_Hash,_Pred,_Alloc> type;
#endif
  };
};

#endif
