// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_CONVOLUTION__HPP__
#define __CICADA__NEURON_CONVOLUTION__HPP__ 1

#include <cicada/neuron/layer.hpp>

namespace cicada
{
  namespace neuron
  {
    class Convolution : public Layer
    {
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      
    public:
      tensor_type weight;
      tensor_type bias;
      
      size_type frame_input;
      size_type frame_output;
      size_type kW;
      size_type dW;
    };
  };
};

#endif
