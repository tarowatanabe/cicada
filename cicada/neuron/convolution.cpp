//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>

#include "cicada/neuron/convolution.hpp"

namespace cicada
{
  namespace neuron
  {
    Convolution::Convolution(size_type __frame, size_type __kW, size_type __dW)
      : frame(__frame), kW(__kW), dW(__dW)
    {
      weight = tensor_type(frame, 1);
      bias   = tensor_type(frame, 1);
    }
    
    Convolution::Convolution(const tensor_type& __weight, const tensor_type& __bias, size_type __kW, size_type __dW)
      : weight(__weight), bias(__bias), kW(__kW), dW(__dW)
    {
      frame = weight.rows();

      if (weight.rows() != bias.rows())
	throw std::runtime_error("invalid weight/bias");
      if (weight.cols() != bias.cols())
	throw std::runtime_error("invalid weight/bias");
      if (weight.cols() != 1)
	throw std::runtime_error("invalid weight");
      if (bias.cols() != 1)
	throw std::runtime_error("invalid bias");
    }

    void Convolution::forward(const tensor_type& data_input)
    {
      const size_type n_frame_input  = data_input.cols();
      const size_type n_frame_output = (n_frame_input - kW) / dW + 1;
      
      data_output.resize(frame, n_frame_output);
      
      for (size_type k = 0; k != n_frame_output; ++ k)
	data_output.col(k) = data_input.block(0, k * dW, data_input.rows(), kW).rowwise().sum().array() * weight.array() + bias.array();
    }
    
    void Convolution::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.resizeLike(data_input);
      gradient_input.setZero();
      
      const size_type n_frame_output = gradient_output.cols();
      
      for (size_type k = 0; k != n_frame_output; ++ k)
	gradient_input.block(0, k * dW, data_input.rows(), kW) += (gradient_output.col(k).transpose() * weight).replicate(1, kW);
    }

    void Convolution::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      const size_type n_frame_output = gradient_output.cols();

      gradient_weight.resizeLike(weight);
      gradient_bias.resizeLike(bias);

      gradient_weight.setZero();
      gradient_bias.setZero();
      
      for (size_type k = 0; k != n_frame_output; ++ k) {
	gradient_weight.array() += data_input.block(0, k * dW, data_input.rows(), kW).rowwise().sum().array() * gradient_output.col(k).array();
	gradient_bias   += gradient_output.col(k);
      }
    }
  }
}
