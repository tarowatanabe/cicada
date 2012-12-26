//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>
#include <cmath>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/tanh.hpp"

namespace cicada
{
  namespace neuron
  {
    float tanh(float x)
    {
      return std::tanh(x);
    }

    void Tanh::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().unaryExpr(std::ptr_fun(tanh));
    }
    
    void Tanh::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = - gradient_output.array() * (data_output.array() * data_output.array() - 1.0);
    }

    std::ostream& Tanh::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os), karma::lit("{\"neuron\":\"tanh\"}"));
      
      return os;
    }
  }
}
