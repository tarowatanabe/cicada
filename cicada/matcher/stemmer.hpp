// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MATCHER__STEMMER__HPP__
#define __CICADA__MATCHER__STEMMER__HPP__ 1

#include <cicada/stemmer.hpp>
#include <cicada/matcher.hpp>

namespace cicada
{
  namespace matcher
  {
    class Stemmer : public cicada::Matcher
    {
    public:
      Stemmer(const cicada::Stemmer* __stemmer) : stemmer(__stemmer) {}
      
    public:
      bool operator()(const symbol_type& x, const symbol_type& y) const
      {
	return stemmer->operator()(x) == stemmer->operator()(y);
      }
      
    private:
      const cicada::Stemmer* stemmer;
    };
  };
};

#endif
