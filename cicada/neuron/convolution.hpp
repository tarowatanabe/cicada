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
      Convolution(size_type frame, size_type kW=1, size_type dW=1);
      Convolution(const tensor_type& weight, const tensor_type& bias, size_type kW=1, size_type dW=1);
      
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual void accumulate(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual layer_ptr_type clone() const { return layer_ptr_type(new Convolution(*this)); }
      
    public:
      tensor_type weight;
      tensor_type bias;

      tensor_type gradient_weight;
      tensor_type gradient_bias;
      
      size_type frame;
      size_type kW;
      size_type dW;
    };
  };
};

#endif
