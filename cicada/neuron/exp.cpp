//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/exp.hpp"

namespace cicada
{
  namespace neuron
  {
    void Exp::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().exp();
    }
    
    void Exp::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.array() = gradient_output.array() * data_output.array();
    }

    std::ostream& Exp::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os), karma::lit("{\"neuron\":\"exp\"}"));
      
      return os;
    }
  }
}
