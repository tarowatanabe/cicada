// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__HEAD__COLLINS__HPP__
#define __CICADA__HEAD__COLLINS__HPP__ 1

namespace cicada
{
  namespace head
  {
    class Collins : public cicada::HeadFinder
    {
    public:
      Collins() : HeadFinder("collins") {}
      
      size_type find_marked_head(const rule_type& rule, const symbol_type& parent) const { return size_type(-1); }
      size_type find_head(const rule_type& rule, const symbol_type& parent) const
      {
	
	
      }
    };
  };
};

#endif
