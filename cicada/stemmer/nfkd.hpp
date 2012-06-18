// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_NFKD__HPP__
#define __CICADA__STEMMER_NFKD__HPP__ 1

#include <vector>

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class NFKD : public Stemmer
    {
    private:
      typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;

    public:
      NFKD();
      
    public:
      symbol_type operator[](const symbol_type& x) const;
      
    private:
      symbol_set_type cache;
      const void* handle;
    };
  };
};

#endif
