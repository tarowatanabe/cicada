//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/min.hpp"

namespace cicada
{
  namespace neuron
  {
    void Min::forward(const tensor_type& data_input)
    {
      if (dimension) {
	// select row and compute column-min
	data_output.conservativeResize(data_input.rows(), 1);
	indices.resize(data_input.rows());
	
	for (size_type row = 0; row != data_input.rows(); ++ row) {
	  int col_min = 0;
	  data_output.col(0)[row] = data_input.row(row).minCoeff(&col_min);
	  indices[row] = col_min;
	}
      } else {
	// select column and compute row-min!
	data_output.conservativeResize(data_input.cols(), 1);
	indices.resize(data_input.cols());
	
	for (size_type col = 0; col != data_input.cols(); ++ col) {
	  int row_min = 0;
	  data_output.col(0)[col] = data_input.col(col).minCoeff(&row_min);
	  indices[col] = row_min;
	}
      }
    }
    
    void Min::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.conservativeResizeLike(data_input);
      gradient_input.setZero();
      
      if (dimension) {
	for (int row = 0; row != data_input.rows(); ++ row)
	  gradient_input.row(row)[indices[row]] = gradient_output.col(0)[row];
      } else {
	for (int col = 0; col != data_input.cols(); ++ col)
	  gradient_input.col(col)[indices[col]] = gradient_output.col(0)[col];
      }
    }
  }
}
