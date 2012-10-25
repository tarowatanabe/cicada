// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_ARABIC__HPP__
#define __CICADA__STEMMER_ARABIC__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    struct ArabicImpl;

    class Arabic : public Stemmer
    {
    private:
      typedef ArabicImpl impl_type;
      
    public:
      Arabic();
      ~Arabic();
      
    public:
      std::string operator()(const utils::piece& word) const;
      
    private:
      impl_type* pimpl;
    };
  };
};

#endif
