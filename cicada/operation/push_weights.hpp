// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__PUSH_WEIGHTS__HPP__
#define __CICADA__OPERATION__PUSH_WEIGHTS__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {

    class PushWeights : public cicada::Operation
    {
    public:
      PushWeights(const std::string& parameter, const int __debug);
  
      void operator()(data_type& data) const;

      bool left;
      bool frontier;
      bool root;
  
      int debug;
    };

  };
};

#endif
