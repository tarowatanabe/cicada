// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_ARABIC__HPP__
#define __CICADA__STEMMER_ARABIC__HPP__ 1

#include <vector>

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    struct ArabicImpl;

    class Arabic : public Stemmer
    {
    private:
      typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
      typedef ArabicImpl impl_type;
      
    public:
      Arabic();
      ~Arabic();
      
    public:
      symbol_type operator[](const symbol_type& x) const;
      
    private:
      impl_type* pimpl;
      symbol_set_type cache;
    };
  };
};

#endif
