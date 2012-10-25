// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_LATIN__HPP__
#define __CICADA__STEMMER_LATIN__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    struct LatinImpl;
    
    class Latin : public Stemmer
    {
    private:
      typedef LatinImpl impl_type;
    
    public:
      Latin();
      ~Latin();
    
    public:
      std::string operator()(const utils::piece& word) const;
    
    private:
      impl_type* pimpl;
    };
  };
};

#endif
