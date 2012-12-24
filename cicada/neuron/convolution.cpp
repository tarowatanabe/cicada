//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/convolution.hpp"

namespace cicada
{
  namespace neuron
  {
    void Convolution::forward(const tensor_type& data_input)
    {
      const size_type n_frame_input  = data_input.cols();
      const size_type n_frame_output = (n_frame_input - kW) / dW + 1;
      
      data_output.resize(frame_ouput, n_frame_output);
      
    }
    
    void Convolution::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      
    }
  }
}
