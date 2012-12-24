//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/softmax.hpp"

namespace cicada
{
  namespace neuron
  {
    void SoftMax::forward(const tensor_type& data_input)
    {
      // use of redux for computing logsum...
      
      const double infty = - std::numeric_limits<double>::infinity();
      
      double logsum = infty;
      for (size_type i = 0; i != data_input.size(); ++ i) {
	const double value = data_input[i];
	
	if (logsum == infty)
	  logsum = value;
	else if (value > infty) {
	  if (logsum >= value)
	    logsum = logsum + utils::mathop::log1p(std::exp(value - logsum));
	  else
	    logsum = value  + utils::mathop::log1p(std::exp(logsum - value));
	}
      }
      
      data_output = (data_input.array() - logsum).exp();
    }
    
    void SoftMax::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      const double sum = (gradient_output.array() * data_output.array()).sum();
      
      gradient_input = (gradient_output.array() - sum) * data_output.array();
    }
  }
}
