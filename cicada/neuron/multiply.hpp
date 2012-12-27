// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_MULTIPLY__HPP__
#define __CICADA__NEURON_MULTIPLY__HPP__ 1

#include <cmath>

#include <cicada/neuron/layer.hpp>

namespace cicada
{
  namespace neuron
  {
    class Multiply : public Layer
    {
    public:
      Multiply();
      Multiply(const tensor_type& weight);
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual void accumulate(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual layer_ptr_type clone() const { return layer_ptr_type(new Multiply(*this)); }
      virtual std::ostream& write(std::ostream& os) const;
      
    public:
      tensor_type weight;
      tensor_type gradient_weight;
    };
  };
};

#endif
