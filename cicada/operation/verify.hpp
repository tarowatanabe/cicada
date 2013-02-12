// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__VERIFY__HPP__
#define __CICADA__OPERATION__VERIFY__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {

    class Verify : public cicada::Operation
    {
    public:
      Verify(const std::string& parameter, const int __debug);
  
      void operator()(data_type& data) const;
      
      int debug;
    };

  };
};

#endif
