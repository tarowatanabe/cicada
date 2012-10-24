// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_KATAKANA__HPP__
#define __CICADA__STEMMER_KATAKANA__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    struct KatakanaImpl;
    
    class Katakana : public Stemmer
    {
    private:
      typedef KatakanaImpl impl_type;
    
    public:
      Katakana();
      ~Katakana();
    
    public:
      std::string operator()(const utils::piece& word) const;
    
    private:
      impl_type*           pimpl;
    };
  };
};

#endif
