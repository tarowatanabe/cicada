// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_SNOWBALL__HPP__
#define __CICADA__STEMMER_SNOWBALL__HPP__ 1

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    struct SnowballImpl;

    class Snowball : public Stemmer
    {
    private:
      typedef SnowballImpl impl_type;
      
    public:
      Snowball(const std::string& language);
      ~Snowball();
    
    public:
      std::string operator()(const utils::piece& word) const;
    
    private:
      impl_type* pimpl;
    };
  };
};

#endif
