// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SIGNATURE_ENGLISH__HPP__
#define __CICADA__SIGNATURE_ENGLISH__HPP__ 1

#include <vector>

#include <cicada/signature.hpp>

namespace cicada
{
  namespace signature
  {
    class English : public Signature
    {
    private:
      typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
      
    public:
      English();
      ~English();
      
    public:
      symbol_type operator[](const symbol_type& x) const;
      
    private:
      symbol_set_type cache;
      void*           pimpl;
    };
  };
};

#endif
