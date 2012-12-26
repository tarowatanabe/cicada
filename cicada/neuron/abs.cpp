//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/abs.hpp"

namespace cicada
{
  namespace neuron
  {
    float abs_derivative(float x)
    {
      return (x >= 0.0 ? 1.0 : -1.0);
    }
    
    void Abs::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().abs();
    }
    
    void Abs::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = (data_input.array().unaryExpr(std::ptr_fun(abs_derivative))) * gradient_output.array();
    }

    std::ostream& Abs::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os), karma::lit("{\"neuron\":\"abs\"}"));
      
      return os;
    }
  }
}
