// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_LOWER__HPP__
#define __CICADA__STEMMER_LOWER__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class Lower : public Stemmer
    {
    public:
      Lower();
      ~Lower();
      
    public:
      std::string operator()(const utils::piece& word) const;
      
    private:
      void* pimpl;
    };
  };
};

#endif
