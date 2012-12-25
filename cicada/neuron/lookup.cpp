//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/lookup.hpp"

namespace cicada
{
  namespace neuron
  {
    void Lookup::forward(const tensor_type& data_input)
    {
      data_output.resize(size, data_input.rows());
      
      // we will cast data_input as a sequence of index...
      const uint32_t* first = reinterpret_cast<const uint32_t*>(data_input.data());
      const uitn32_t* last = first + data_input.rows();
      for (size_type frame = 0; first != last; ++ first, ++ frame) {
	// perform resizing of weight
	if (*first >= weight.cols()) {
	  // from weight.cols() to *first + 1, assign random value...!
	  const size_type frame_first = weight.cols();
	  const size_type frame_last  = *first + 1;
	  
	  weight.conservativeResize(Eigen::NoChange_t, *first + 1);
	  for (size_type i = frame_first; i != frame_last; ++ i)
	    weight.cols(i).setRandom();
	}
	
	data_output.col(frame) = weight.col(*first);
      }
    }
    
    void Lookup::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      
    }
    
    void Lookup::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      
    }
  }
}
