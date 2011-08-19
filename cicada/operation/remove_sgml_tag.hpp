// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__REMOVE_SGML_TAG__HPP__
#define __CICADA__OPERATION__REMOVE_SGML_TAG__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {

    class RemoveSGMLTag : public cicada::Operation
    {
    public:
      RemoveSGMLTag(const std::string& parameter, const int __debug);
  
      void operator()(data_type& data) const;
  
      bool lattice_mode;
      bool forest_mode;
      bool remove_bos_eos;
      int debug;
    };

  };
};

#endif
