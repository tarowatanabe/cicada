//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/copy.hpp"

namespace cicada
{
  namespace neuron
  {
    void Copy::forward(const tensor_type& data_input)
    {
      data_output = data_input;
    }
    
    void Copy::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = gradient_output;
    }

    std::ostream& Copy::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os), karma::lit("{\"neuron\":\"copy\"}"));
      
      return os;
    }
  }
}
