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
      
      data_output.resize(frame_output, n_frame_output);
      
      for (size_type k = 0; k != n_frame_output; ++ k)
	data_output.col(k) = data_input.block(, ).colwise().sum() * weight + bias;
    }
    
    void Convolution::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.resizeLike(data_input);
      
      const size_type n_frame_output = gradient_output.cols();
      
      for (size_type k = 0; k != n_frame_output; ++ k) {
	
      }
    }
  }
}
