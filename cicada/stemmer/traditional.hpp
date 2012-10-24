// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_TRADITIONAL__HPP__
#define __CICADA__STEMMER_TRADITIONAL__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class Traditional : public Stemmer
    {
    public:
      Traditional();
      ~Traditional();
      
    public:
      std::string operator()(const utils::piece& word) const;
      
    private:
      void* pimpl;
    };
  };
};

#endif
