// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__CLEAR__HPP__
#define __CICADA__OPERATION__CLEAR__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {
    class Clear : public Operation
    {
    public:
      Clear(const std::string& parameter,
	    const int __debug);
  
      void operator()(data_type& data) const;
      
      bool clear_hypergraph;
      bool clear_lattice;
      bool clear_spans;
      bool clear_targets;
      bool clear_counts;
      
      int debug;
    };

  };
};


#endif
