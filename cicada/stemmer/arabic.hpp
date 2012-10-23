// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_ARABIC__HPP__
#define __CICADA__STEMMER_ARABIC__HPP__ 1

#include <cicada/stemmer.hpp>

#include <utils/array_power2.hpp>

namespace cicada
{
  namespace stemmer
  {
    struct ArabicImpl;

    class Arabic : public Stemmer
    {
    private:
      typedef std::pair<symbol_type, symbol_type> symbol_pair_type;
      typedef utils::array_power2<symbol_pair_type, 1024 * 8, std::allocator<symbol_pair_type> > symbol_pair_set_type;
      typedef ArabicImpl impl_type;
      
    public:
      Arabic();
      ~Arabic();
      
    public:
      symbol_type operator[](const symbol_type& x) const;
      
    private:
      impl_type* pimpl;
      symbol_pair_set_type cache;
    };
  };
};

#endif
