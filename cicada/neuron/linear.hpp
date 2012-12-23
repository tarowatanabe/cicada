// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_LINEAR__HPP__
#define __CICADA__NEURON_LINEAR__HPP__ 1

#include <cicada/neuron/layer.hpp>

namespace cicada
{
  namespace neuron
  {
    class Linear : public Layer
    {
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      
    public:
      tensor_type weight;
      tensor_type bias;
    };
  };
};

#endif
