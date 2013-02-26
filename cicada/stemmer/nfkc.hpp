// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
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
      ~NFKC();
      
    public:
      std::string operator()(const utils::piece& word) const;
      
    private:
      void* handle;
    };
  };
};

#endif
