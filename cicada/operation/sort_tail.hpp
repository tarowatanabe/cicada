// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__SORT_TAIL__HPP__
#define __CICADA__OPERATION__SORT_TAIL__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {
    class SortTail : public Operation
    {
    public:
      SortTail(const std::string& parameter, const int __debug);

      void operator()(data_type& data) const;
      
      int debug;
    };

  };
};


#endif
