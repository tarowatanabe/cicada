// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MATCHER__LOWER__HPP__
#define __CICADA__MATCHER__LOWER__HPP__ 1

#include <cicada/stemmer.hpp>
#include <cicada/matcher.hpp>

namespace cicada
{
  namespace matcher
  {
    class Lower : public cicada::Matcher
    {
    public:
      Lower() : lower(&cicada::Stemmer::create("lower")) {}
      
    public:
      bool operator()(const symbol_type& x, const symbol_type& y) const
      {
	return lower->operator()(x) == lower->operator()(y);
      }
      
    private:
      const cicada::Stemmer* lower;
    };
  };
};

#endif
