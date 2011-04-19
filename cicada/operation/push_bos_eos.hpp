// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__PUSH_BOS_EOS__HPP__
#define __CICADA__OPERATION__PUSH_BOS_EOS__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {

    class PushBosEos : public cicada::Operation
    {
    public:
      PushBosEos(const std::string& parameter, const int __debug);
  
      void operator()(data_type& data) const;
  
      int debug;
    };

  };
};

#endif
