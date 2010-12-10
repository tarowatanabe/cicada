// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MATCHER__WORDNET__HPP__
#define __CICADA__MATCHER__WORDNET__HPP__ 1

#include <string>

#include <cicada/matcher.hpp>

namespace cicada
{
  namespace matcher
  {
    class WordNet : public cicada::Matcher
    {
    public:
      WordNet();
      WordNet(const std::string& path);
      
    public:
      bool operator()(const symbol_type& x, const symbol_type& y) const;
    };
  };
};

#endif
