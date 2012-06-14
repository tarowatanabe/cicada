// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__STREE_MAP__HPP__
#define __UTILS__STREE_MAP__HPP__ 1

#include <utils/smap.hpp>

namespace utils
{
   template <typename _Key,
	     typename _Tp,
	     typename _Pred=std::less<_Key>,
	     typename _Alloc=std::allocator<std::pair<const _Key, _Tp> >,
	     size_t BN=64>
   struct stree_map
   {
     typedef sti::smap<_Key,_Tp,_Pred, _Aloc, BN> type;
   };
};

#endif
