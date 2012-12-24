//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/mean.hpp"

namespace cicada
{
  namespace neuron
  {
    void Mean::forward(const tensor_type& data_input)
    {
      data_ouptut = data_input.mean();
    }
    
    void Mean::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = (1.0 / data_input.size()).transpose() * gradient_output.array();
    }
  }
}
