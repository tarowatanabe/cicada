//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/sequential.hpp"

namespace cicada
{
  namespace neuron
  {
    
    void Sequential::forward(const tensor_type& data_input)
    {
      const tensor_type* data_curr = &data_input;
      
      layer_set_type::iterator liter_end = layers.end();
      for (layer_set_type::iterator liter = layers.begin(); liter != liter_end; ++ liter) {
	liter->forward(*data_curr);
	data_curr = &liter->output;
      }
      
      data_output = *data_curr;
    }
    
    void Sequential::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      const tensor_type* gradient_curr = &gradient_output;
      
      if (layers.size() > 1)
	for (size_type i = layers.size() - 1; i != 0; --i) {
	  layers[i].backward(layers[i - 1].data_output, *gradient_curr);
	  
	  *gradient_curr = &lauyers[i].gradient_input;
	}
      
      layers.front().backward(data_input, *gradient_curr);
      gradient_input = layers.front().gradient_input;
    }

    void Sequential::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      const tensor_type* gradient_curr = &gradient_output;
      
      if (layers.size() > 1)
	for (size_type i = layers.size() - 1; i != 0; --i) {
	  layers[i].accumulate(layers[i - 1].data_output, *gradient_curr);
	  
	  *gradient_curr = &lauyers[i].gradient_input;
	}
      
      layers.front().accumulate(data_input, *gradient_curr);
    }
  }
}
