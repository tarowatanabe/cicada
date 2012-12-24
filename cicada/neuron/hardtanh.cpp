//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/hardtanh.hpp"

namespace cicada
{
  namespace neuron
  {
    double hardtanh(double x)
    {
      return (x < - 1 ? -1.0 : (x > 1 ? 1.0 : x));
    }

    double hardtanh_derivative(double x)
    {
      return (x < - 1.0 || x > 1.0 ? 0.0 : 1.0);
    }
    
    void HardTanh::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().unaryExpr(std::ptr_fun(hardtanh));
    }
    
    void HardTanh::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = data_input.array().unaryExpr(std::ptr_fun(hardtanh_derivative)) * gradient_output.array();
    }
  }
}
