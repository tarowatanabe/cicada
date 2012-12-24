//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/log.hpp"

namespace cicada
{
  namespace neuron
  {
    void Log::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().log();
    }
    
    void Log::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.array() = data_input.array() / gradient_output.array();
    }
  }
}
