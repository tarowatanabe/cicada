//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/exp.hpp"

namespace cicada
{
  namespace neuron
  {
    void Exp::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().exp();
    }
    
    void Exp::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.array() = gradient_output.array() * data_output.array();
    }
  }
}
