//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <cmath>

#include "cicada/neuron/tanh.hpp"

namespace cicada
{
  namespace neuron
  {
    float tanh(float x)
    {
      return std::tanh(x);
    }

    void Tanh::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().unaryExpr(std::ptr_fun(tanh));
    }
    
    void Tanh::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = - gradient_output.array() * (data_output.array() * data_output.array() - 1.0);
    }
  }
}
