//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/square.hpp"

namespace cicada
{
  namespace neuron
  {
    void Square::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().square();
    }
    
    void Square::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = 2.0 * gradient_output.array() * data_input.array();
    }

    std::ostream& Square::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os), karma::lit("{\"neuron\":\"square\"}"));
      
      return os;
    }
  }
}
