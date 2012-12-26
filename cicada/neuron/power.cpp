//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/power.hpp"

namespace cicada
{
  namespace neuron
  {
    void Power::forward(const tensor_type& data_input)
    {
      data_output = data_input.array().pow(pow);
    }
    
    void Power::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = pow * gradient_output.array() * data_output.array() / data_input.array();
    }
    
    struct power_policy : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double)
      {
	return 20;
      }
    };

    std::ostream& Power::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;

      karma::real_generator<double, power_policy> double20;
      
      karma::generate(std::ostream_iterator<char>(os),
		      karma::lit('{')
		      << karma::lit("\"neuron\"") << ':' << karma::lit("\"power\"")
		      << ',' << "\"power\"" << ':' << double20
		      << '}',
		      pow);
      
      return os;
    }
  }
}
