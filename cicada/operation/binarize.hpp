// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__BINARIZE__HPP__
#define __CICADA__OPERATION__BINARIZE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {

    class Binarize : public cicada::Operation
    {
    public:
      Binarize(const std::string& parameter, const int __debug);
  
      void operator()(data_type& data) const;
  
      int order;
  
      bool left;
      bool right;
      bool all;
  
      int debug;
    };

  };
};

#endif
