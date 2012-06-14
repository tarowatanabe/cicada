// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__STREE_SET__HPP__
#define __UTILS__STREE_SET__HPP__ 1

#include <utils/sset.hpp>

namespace utils
{
   template <typename _Key,
	     typename _Pred=std::less<_Key>,
	     typename _Alloc=std::allocator<_Key >,
	     size_t BN=64>
   struct stree_set
   {
     typedef sti::sset<_Key,_Pred, _Aloc, BN> type;
   };
};

#endif
