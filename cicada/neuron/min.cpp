//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/min.hpp"

namespace cicada
{
  namespace neuron
  {
    void Min::forward(const tensor_type& data_input)
    {
      data_ouptut = data_input.minCoeff(&index);
    }
    
    void Min::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      // index-pos is one, others zero...
      gradient_input = (1.0).transpose() * gradient_output;
    }
  }
}
