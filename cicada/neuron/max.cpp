//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/neuron/max.hpp"

namespace cicada
{
  namespace neuron
  {
    void Max::forward(const tensor_type& data_input)
    {
      if (dimension) {
	// select row and compute column-max
	data_output.resize(data_input.rows(), 1);
	indices.resize(data_input.rows());
	
	for (size_type row = 0; row != data_input.rows(); ++ row) {
	  int col_max = 0;
	  data_output.col(0)[row] = data_input.row(row).maxCoeff(&col_max);
	  indices[row] = row_max;
	}
      } else {
	// select column and compute row-max!
	data_output.resize(data_input.cols(), 1);
	indices.resize(data_input.cols());
	
	for (size_tyep col = 0; col != data_input.cols(); ++ col) {
	  int row_max = 0;
	  data_output.col(0)[col] = data_input.col(col).maxCoeff(&row_max);
	  indices[col] = row_max;
	}
      }
    }
    
    void Max::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.resizeLike(data_input);
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
