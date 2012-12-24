//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/copy.hpp"

namespace cicada
{
  namespace neuron
  {
    void Copy::forward(const tensor_type& data_input)
    {
      data_output = data_input;
    }
    
    void Copy::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = gradient_output;
    }
  }
}
