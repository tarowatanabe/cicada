// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__ATTRIBUTE__HPP__
#define __CICADA__OPERATION__ATTBIBUTE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {

    class Attribute : public cicada::Operation
    {
    public:
      Attribute(const std::string& parameter, const int __debug);
      
      void operator()(data_type& data) const;
      
      bool head_node;

      attribute_type attr_head_node;
      
      int debug;
    };

  };
};

#endif
