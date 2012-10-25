// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_HALFWIDTH__HPP__
#define __CICADA__STEMMER_HALFWIDTH__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class Halfwidth : public Stemmer
    {
    public:
      Halfwidth();
      ~Halfwidth();
      
    public:
      std::string operator()(const utils::piece& word) const;
      
    private:
      void* pimpl;
    };
  };
};

#endif
