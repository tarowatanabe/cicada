// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_SUFFIX__HPP__
#define __CICADA__STEMMER_SUFFIX__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class Suffix : public Stemmer
    {
    public:
      Suffix(const size_type __size) : size(__size) {}
    
    public:
      std::string operator()(const utils::piece& word) const;    
      
    private:
      size_type size;
    };
  };
};

#endif
