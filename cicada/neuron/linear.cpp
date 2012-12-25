//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>

#include "cicada/neuron/linear.hpp"

namespace cicada
{
  namespace neuron
  {
    Linear::Linear(size_type size_input, size_type size_output)
    {
      if (! size_input)
	throw std::runtime_error("invalid input size");
      if (! size_output)
	throw std::runtime_error("invalid output size");
      
      weight = tensor_type::Random(size_output, size_input);
      bias   = tensor_type::Random(size_output, 1);
    }
    
    Linear::Linear(const tensor_type& __weight, const tensor_type& __bias)
      : weight(__weight), bias(__bias)
    {
      if (bias.cols() != 1)
	throw std::runtime_error("invalid bias");
      
      if (bias.rows() != weight.rows())
	throw std::runtime_error("invalid weight");
    }

    void Linear::forward(const tensor_type& data_input)
    {
      if (data_input.rows() != weight.cols())
	throw std::runtime_error("invalid input");
      
      data_output = bias + weight * data_input;
    }
    
    void Linear::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      if (data_input.rows() != weight.cols())
	throw std::runtime_error("invalid input");

      if (weight.rows() != gradient_output.rows())
	throw std::runtime_error("invalid gradient output");
      
      gradient_input = (weight.transpose() * gradient_output).array() + 1.0;
    }
    
    void Linear::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      if (data_input.rows() != weight.cols())
	throw std::runtime_error("invalid input");
      
      if (weight.rows() != gradient_output.rows())
	throw std::runtime_error("invalid gradient output");
      
      gradient_weight = gradient_output * data_input.transpose();
      gradient_bias   = gradient_output;
    }
  }
}
