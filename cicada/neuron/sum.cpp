//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/sum.hpp"

namespace cicada
{
  namespace neuron
  {
    void Sum::forward(const tensor_type& data_input)
    {
      if (dimension) {
	// select row and compute column-sum
	data_output.resize(data_input.rows(), 1);
	
	for (size_type row = 0; row != data_input.rows(); ++ row)
	  data_output.col(0)[row] = data_input.row(row).sum();
      } else {
	// select column and compute row-sum
	data_output.resize(data_input.cols(), 1);
	
	for (size_type col = 0; col != data_input.cols(); ++ col)
	  data_output.col(0)[col] = data_input.col(col).sum();
      }
    }
    
    void Sum::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.resizeLike(data_input);
      
      if (dimension) {
	for (int row = 0; row != data_input.rows(); ++ row)
	  gradient_input.row(row).setConstant(gradient_output.col(0)[row]);
      } else {
	for (int col = 0; col != data_input.cols(); ++ col)
	  gradient_input.col(col).setConstant(gradient_output.col(0)[col]);
      }
    }
  }
}
