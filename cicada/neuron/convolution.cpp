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
      
      data_output.conservativeResize(frame, n_frame_output);
      
      for (size_type k = 0; k != n_frame_output; ++ k)
	data_output.col(k) = data_input.block(0, k * dW, data_input.rows(), kW).rowwise().sum().array() * weight.array() + bias.array();
    }
    
    void Convolution::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.conservativeResizeLike(data_input);
      gradient_input.setZero();
      
      const size_type n_frame_output = gradient_output.cols();
      
      for (size_type k = 0; k != n_frame_output; ++ k)
	gradient_input.block(0, k * dW, data_input.rows(), kW) += (gradient_output.col(k).transpose() * weight).replicate(1, kW);
    }

    void Convolution::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      const size_type n_frame_output = gradient_output.cols();

      gradient_weight.conservativeResizeLike(weight);
      gradient_bias.conservativeResizeLike(bias);

      gradient_weight.setZero();
      gradient_bias.setZero();
      
      for (size_type k = 0; k != n_frame_output; ++ k) {
	gradient_weight.array() += data_input.block(0, k * dW, data_input.rows(), kW).rowwise().sum().array() * gradient_output.col(k).array();
	gradient_bias   += gradient_output.col(k);
      }
    }
  }
}
