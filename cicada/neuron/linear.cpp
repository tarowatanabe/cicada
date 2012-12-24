

#include "cicada/neuron/linear.hpp"

namespace cicada
{
  namespace neuron
  {
    void Linear::forward(const tensor_type& data_input)
    {
      // assume eigen's operator!
      data_output = bias + weight * data_input;
    }
    
    void Linear::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      // assume eigen's operator!
      gradient_input = weight.transpose() * gradient_output + 1;
    }
  }
}
