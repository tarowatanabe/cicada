//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/concat.hpp"

#include <stdexcept>
#include <memory>

namespace cicada
{
  namespace neuron
  {
    void Concat::forward(const tensor_type& data_input)
    {
      sizes.resize(layers.size());
      data_output.resize(0, 0);

      if (dimension) {
	size_set_type::iterator siter = sizes.begin();
	layer_set_type::iterator liter_end = layers.end();
	for (layer_set_type::iterator liter = layers.begin(); liter != liter_end; ++ liter, ++ siter) {
	  (*liter)->forward(data_input);
	  
	  *siter = (*liter)->data_output.rows();
	  
	  if (data_output.cols()) {
	    if (data_output.rows() != static_cast<int>(*siter))
	      throw std::runtime_error("invalid concat");

	    data_output.conservativeResize(Eigen::NoChange, data_output.cols() + 1);
	    data_output.col(data_output.cols() - 1) = data_output.col(0);
	  } else
	    data_output = (*liter)->data_output.col(0);
	}
      } else {
	size_set_type::iterator siter = sizes.begin();
	layer_set_type::iterator liter_end = layers.end();
	for (layer_set_type::iterator liter = layers.begin(); liter != liter_end; ++ liter, ++ siter) {
	  (*liter)->forward(data_input);
	  
	  *siter = (*liter)->data_output.rows();
	  
	  if (data_output.rows()) {
	    data_output.conservativeResize(data_output.rows() + *siter, 1);
	    data_output.block(data_output.rows() - *siter, 0, *siter, 1) = data_output.col(0);
	  } else
	    data_output = (*liter)->data_output.col(0);
	}
      }
    }
    
    void Concat::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.conservativeResizeLike(data_input);
      gradient_input.setZero();
      
      if (dimension) {
	for (size_type i = 0; i != sizes.size(); ++ i) {
	  layers[i]->backward(data_input, gradient_output.col(i));
	  gradient_input += layers[i]->gradient_input;
	}
      } else {
	size_type offset = 0;
	for (size_type i = 0; i != sizes.size(); ++ i) {
	  layers[i]->backward(data_input, gradient_output.block(offset, 0, sizes[i], 1));
	  gradient_input += layers[i]->gradient_input;
	  offset += sizes[i];
	}
      }
    }

    void Concat::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      if (dimension) {
	for (size_type i = 0; i != sizes.size(); ++ i)
	  layers[i]->accumulate(data_input, gradient_output.col(i));
      } else {
	size_type offset = 0;
	for (size_type i = 0; i != sizes.size(); ++ i) {
	  layers[i]->accumulate(data_input, gradient_output.block(offset, 0, sizes[i], 1));
	  offset += sizes[i];
	}
      }
    }

    Concat::layer_ptr_type Concat::clone() const
    {
      std::auto_ptr<Concat> cloned(new Concat(*this));
      
      for (size_type i = 0; i != layers.size(); ++ i)
	cloned->layers[i] = layers[i]->clone();
      
      return layer_ptr_type(cloned.release());
    }
  }
}
