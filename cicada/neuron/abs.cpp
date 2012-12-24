//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/abs.hpp"

namespace cicada
{
  namespace neuron
  {
    double abs_derivative(double x)
    {
      return (x >= 0.0 ? 1.0 : -1.0);
    }
    
    void Abs::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().abs();
    }
    
    void Abs::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.array() = (data_input.array().unaryExpr(std::ptr_fun(abs_derivative))) * gradient_output.array();
    }
  }
}
