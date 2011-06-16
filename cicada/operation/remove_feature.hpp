// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__REMOVE_FEATURE__HPP__
#define __CICADA__OPERATION__REMOVE_FEATURE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {

    class RemoveFeature : public cicada::Operation
    {
    public:
      RemoveFeature(const std::string& parameter, const int __debug);
      
      void operator()(data_type& data) const;

      typedef feature_function_type::feature_type feature_type;
      
      typedef std::vector<feature_type, std::allocator<feature_type> > remove_set_type;
      
      remove_set_type removes;

      int debug;
    };

  };
};

#endif
