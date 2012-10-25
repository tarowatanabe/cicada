// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SIGNATURE_ENGLISH__HPP__
#define __CICADA__SIGNATURE_ENGLISH__HPP__ 1

#include <cicada/signature.hpp>

namespace cicada
{
  namespace signature
  {
    class English : public Signature
    {
    public:
      English();
      ~English();
      
    public:
      virtual std::string operator()(const utils::piece& word) const;
      
    private:
      void* pimpl;
    };
  };
};

#endif
