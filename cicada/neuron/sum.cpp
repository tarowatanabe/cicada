//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/sum.hpp"

namespace cicada
{
  namespace neuron
  {
    void Sum::forward(const tensor_type& data_input)
    {
      if (dimension)
	data_output = data_input.rowwise().sum();
      else
	data_output = data_input.colwise().sum().transpose();
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

    std::ostream& Sum::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os),
		      karma::lit('{')
		      << karma::lit("\"neuron\"") << ':' << karma::lit("\"sum\"")
		      << ',' << "\"dimension\"" << ':' << karma::bool_
		      << '}',
		      dimension);
      
      return os;
    }
  }
}
