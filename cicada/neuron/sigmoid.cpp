//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/sigmoid.hpp"

namespace cicada
{
  namespace neuron
  {
    void Sigmoid::forward(const tensor_type& data_input)
    {
      data_putput = 1.0 / ((- data_input.array()).exp() + 1.0);
    }
    
    void Sigmoid::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = gradient_output.array() * (1.0 - data_output.array()) * data_output.array();
    }
  }
}
