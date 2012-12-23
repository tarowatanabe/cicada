// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_LOOKUP__HPP__
#define __CICADA__NEURON_LOOKUP__HPP__ 1

//
// lookup-table layer which map input, a sequence of id, into a sequence of tensor
//

#include <cicada/neuron/layer.hpp>

namespace cicada
{
  namespace neuron
  {
    class Lookup : public Layer
    {
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      
    public:
      
    };
  };
};

#endif
