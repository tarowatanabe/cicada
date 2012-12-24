//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/power.hpp"

namespace cicada
{
  namespace neuron
  {
    void Power::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().pow(pow);
    }
    
    void Power::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = pow * gradient_output.array() * data_output.array() / data_input.array();
    }
  }
}
