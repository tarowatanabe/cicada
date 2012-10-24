// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_DIGIT__HPP__
#define __CICADA__STEMMER_DIGIT__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class Digit : public Stemmer
    {
    public:
      Digit() {}
    
    public:
      std::string operator()(const utils::piece& word) const;
    };
  };
};

#endif
