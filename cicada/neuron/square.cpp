//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/square.hpp"

namespace cicada
{
  namespace neuron
  {
    void Square::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().square();
    }
    
    void Square::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.array() = 2.0 * gradient_output.array() * data_input.array();
    }
  }
}
