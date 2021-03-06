//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/hardtanh.hpp"

namespace cicada
{
  namespace neuron
  {
    float hardtanh(float x)
    {
      return (x < - 1 ? -1.0 : (x > 1 ? 1.0 : x));
    }
    
    float hardtanh_derivative(float x)
    {
      return (x < - 1.0 || x > 1.0 ? 0.0 : 1.0);
    }
    
    void HardTanh::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().unaryExpr(std::ptr_fun(hardtanh));
    }
    
    void HardTanh::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = data_input.array().unaryExpr(std::ptr_fun(hardtanh_derivative)) * gradient_output.array();
    }

    std::ostream& HardTanh::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os), karma::lit("{\"neuron\":\"hardtanh\"}"));
      
      return os;
    }
  }
}
