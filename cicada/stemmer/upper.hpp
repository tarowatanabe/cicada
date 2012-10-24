// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_UPPER__HPP__
#define __CICADA__STEMMER_UPPER__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class Upper : public Stemmer
    {
    public:
      Upper();
      ~Upper();
      
    public:
      std::string operator()(const utils::piece& word) const;
      
    private:
      void* pimpl;
    };
  };
};

#endif
