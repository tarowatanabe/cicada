// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_HIRAGANA__HPP__
#define __CICADA__STEMMER_HIRAGANA__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    struct HiraganaImpl;
    
    class Hiragana : public Stemmer
    {
    private:
      typedef HiraganaImpl impl_type;
    
    public:
      Hiragana();
      ~Hiragana();
    
    public:
      std::string operator()(const utils::piece& word) const;
    
    private:
      impl_type*           pimpl;
    };
  };
};

#endif
