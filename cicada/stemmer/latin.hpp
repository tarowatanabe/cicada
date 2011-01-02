// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_LATIN__HPP__
#define __CICADA__STEMMER_LATIN__HPP__ 1

#include <vector>

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    struct LatinImpl;
    
    class Latin : public Stemmer
    {
    private:
      typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
      typedef LatinImpl impl_type;
    
    public:
      Latin();
      ~Latin();
    
    public:
      symbol_type operator[](const symbol_type& x) const;
    
    private:
      symbol_set_type cache;
      impl_type*      pimpl;
    };
  };
};

#endif
