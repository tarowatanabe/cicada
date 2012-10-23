// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SIGNATURE_CHINESE__HPP__
#define __CICADA__SIGNATURE_CHINESE__HPP__ 1

#include <cicada/signature.hpp>

#include <utils/array_power2.hpp>

namespace cicada
{
  namespace signature
  {
    class Chinese : public Signature
    {
    private:
      typedef std::pair<symbol_type, symbol_type> symbol_pair_type;
      typedef utils::array_power2<symbol_pair_type, 1024 * 16, std::allocator<symbol_pair_type> > symbol_pair_set_type;
      
    public:
      Chinese();
      ~Chinese();
    public:
      symbol_type operator[](const symbol_type& x) const;
      
    private:
      symbol_pair_set_type cache;
      void* pimpl;
    };
  };
};

#endif
