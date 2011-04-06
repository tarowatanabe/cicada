// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__REMOVE_EPSILON__HPP__
#define __CICADA__OPERATION__REMOVE_EPSILON__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {

    class RemoveEpsilon : public cicada::Operation
    {
    public:
      RemoveEpsilon(const std::string& parameter, const int __debug);
  
      void operator()(data_type& data) const;
  
      bool lattice_mode;
      bool forest_mode;
      int debug;
    };

  };
};

#endif
