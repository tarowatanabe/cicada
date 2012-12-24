//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/sqrt.hpp"

namespace cicada
{
  namespace neuron
  {
    void Sqrt::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().sqrt();
    }
    
    void Sqrt::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = 0.5 * gradient_output.array() / data_output.array();
    }
  }
}
