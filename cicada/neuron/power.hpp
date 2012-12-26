// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_POWER__HPP__
#define __CICADA__NEURON_POWER__HPP__ 1

#include <cmath>

#include <cicada/neuron/layer.hpp>

namespace cicada
{
  namespace neuron
  {
    class Power : public Layer
    {
    public:
      Power(const double __pow=1.0) : pow(__pow) {}
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual void accumulate(const tensor_type& data_input, const tensor_type& gradient_output) {}
      virtual layer_ptr_type clone() const { return layer_ptr_type(new Power(*this)); }
      virtual std::ostream& write(std::ostream& os) const;
    private:
      double pow;
    };
  };
};

#endif
