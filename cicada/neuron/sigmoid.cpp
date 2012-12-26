//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/sigmoid.hpp"

namespace cicada
{
  namespace neuron
  {
    void Sigmoid::forward(const tensor_type& data_input)
    {
      data_output = 1.0 / ((- data_input.array()).exp() + 1.0);
    }
    
    void Sigmoid::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = gradient_output.array() * (1.0 - data_output.array()) * data_output.array();
    }

    std::ostream& Sigmoid::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os), karma::lit("{\"neuron\":\"sigmoid\"}"));
      
      return os;
    }
  }
}
