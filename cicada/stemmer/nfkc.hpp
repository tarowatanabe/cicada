// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_NFKC__HPP__
#define __CICADA__STEMMER_NFKC__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class NFKC : public Stemmer
    {
    public:
      NFKC();
      
    public:
      std::string operator()(const utils::piece& word) const;
      
    private:
      const void* handle;
    };
  };
};

#endif
