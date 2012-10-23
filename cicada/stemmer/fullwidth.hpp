// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STEMMER_FULLWIDTH__HPP__
#define __CICADA__STEMMER_FULLWIDTH__HPP__ 1

#include <cicada/stemmer.hpp>

#include <utils/array_power2.hpp>

namespace cicada
{
  namespace stemmer
  {
    class Fullwidth : public Stemmer
    {
    private:
      typedef std::pair<symbol_type, symbol_type> symbol_pair_type;
      typedef utils::array_power2<symbol_pair_type, 1024 * 8, std::allocator<symbol_pair_type> > symbol_pair_set_type;
      
    public:
      Fullwidth();
      ~Fullwidth();
      
    public:
      symbol_type operator[](const symbol_type& x) const;
      
    private:
      symbol_pair_set_type cache;
      void*                pimpl;
    };
  };
};

#endif
