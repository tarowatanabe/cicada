//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/linear.hpp"

namespace cicada
{
  namespace neuron
  {
    void Linear::forward(const tensor_type& data_input)
    {
      // assume eigen's operator!
      data_output = bias + weight * data_input;
    }
    
    void Linear::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      // assume eigen's operator!
      gradient_input = weight.transpose() * gradient_output + 1;
    }
    
    void Linear::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_weight += gradient_output * data_input.transpose();
      gradient_bias   += gradient_output;
    }
  }
}
