//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/tanh.hpp"

namespace cicada
{
  namespace neuron
  {
    void Tanh::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().unaryExpr(std::ptr_fun(std::tanh));
    }
    
    void Tanh::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.array() = - gradient_output.array() * (data_output.array() * data_output.array() - 1.0);
    }
  }
}
