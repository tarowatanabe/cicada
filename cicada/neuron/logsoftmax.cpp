//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <climits>
#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/logsoftmax.hpp"

#include "utils/mathop.hpp"

namespace cicada
{
  namespace neuron
  {
    void LogSoftMax::forward(const tensor_type& data_input)
    {
      // use of redux for computing logsum...
      
      const double infty = - std::numeric_limits<double>::infinity();
      
      double logsum = infty;
      for (difference_type i = 0; i != data_input.rows(); ++ i) {
	const double value = data_input.col(0)[i];
	
	if (logsum == infty)
	  logsum = value;
	else if (value > infty) {
	  if (logsum >= value)
	    logsum = logsum + utils::mathop::log1p(std::exp(value - logsum));
	  else
	    logsum = value  + utils::mathop::log1p(std::exp(logsum - value));
	}
      }
      
      data_output = data_input.array() - logsum;
    }
    
    void LogSoftMax::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input = gradient_output.array() - (data_output.array().exp() * gradient_output.array().sum());
    }

    std::ostream& LogSoftMax::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os), karma::lit("{\"neuron\":\"log-softmax\"}"));
      
      return os;
    }
  }
}
