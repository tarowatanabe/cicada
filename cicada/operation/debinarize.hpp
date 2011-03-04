// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__DEBINARIZE__HPP__
#define __CICADA__OPERATION__DEBINARIZE__HPP__ 1

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {
    class Debinarize : public Operation
    {
    public:
      Debinarize(const std::string& parameter, const int __debug);
      
      void operator()(data_type& data) const;
  
      int debug;
    };

  };
};


#endif
