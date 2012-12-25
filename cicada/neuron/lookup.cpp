//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>

#include "cicada/neuron/lookup.hpp"

namespace cicada
{
  namespace neuron
  {
    void Lookup::forward(const tensor_type& data_input)
    {
      if (data_input.cols() != 1)
	throw std::runtime_error("invalid input");

      data_output.resize(size, data_input.rows());
      
      // we will cast data_input as a sequence of index...
      const uint32_t* first = reinterpret_cast<const uint32_t*>(data_input.data());
      const uint32_t* last = first + data_input.rows();
      for (size_type frame = 0; first != last; ++ first, ++ frame) {
	// perform resizing of weight
	if (*first >= weight.cols()) {
	  // from weight.cols() to *first + 1, assign random value...!
	  const size_type frame_first = weight.cols();
	  const size_type frame_last  = *first + 1;
	  
	  weight.conservativeResize(size, *first + 1);
	  for (size_type i = frame_first; i != frame_last; ++ i)
	    weight.col(i).setRandom();
	}
	
	data_output.col(frame) = weight.col(*first);
      }
    }
    
    void Lookup::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      // we will do nothing, since this must be the starting point!
      
    }
    
    void Lookup::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      if (data_input.cols() != 1)
	throw std::runtime_error("invalid input");
      
      if (gradient_output.rows() != size)
	throw std::runtime_error("invalid gradient output");
      
      gradient_weight.resizeLike(weight);
      
      // we do not clear all, but only required part...
      //gradient_weight.setZero();
      {
	const uint32_t* first = reinterpret_cast<const uint32_t*>(data_input.data());
	const uint32_t* last = first + data_input.rows();
	for (/**/; first != last; ++ first)
	  gradient_weight.col(*first).setZero();
      }
      
      const uint32_t* first = reinterpret_cast<const uint32_t*>(data_input.data());
      const uint32_t* last = first + data_input.rows();
      for (size_type frame = 0; first != last; ++ first, ++ frame)
	gradient_weight.col(*first) += gradient_output.col(frame);
    }
  }
}
