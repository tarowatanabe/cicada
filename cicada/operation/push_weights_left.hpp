// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__PUSH_WEIGHTS_LEFT__HPP__
#define __CICADA__OPERATION__PUSH_WEIGHTS_LEFT__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {

    class PushWeightsLeft : public cicada::Operation
    {
    public:
      PushWeightsLeft(const std::string& parameter, const int __debug);
  
      void operator()(data_type& data) const;
  
      int debug;
    };

  };
};

#endif
