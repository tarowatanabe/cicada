// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_NFKD__HPP__
#define __CICADA__STEMMER_NFKD__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class NFKD : public Stemmer
    {
    public:
      NFKD();
      
    public:
      std::string operator()(const utils::piece& word) const;
      
    private:
      const void* handle;
    };
  };
};

#endif
