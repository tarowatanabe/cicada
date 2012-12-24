//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/sum.hpp"

namespace cicada
{
  namespace neuron
  {
    void Sum::forward(const tensor_type& data_input)
    {
      data_ouptut = data_input.sum();
    }
    
    void Sum::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = (1.0).transpose() * gradient_output;
    }
  }
}
