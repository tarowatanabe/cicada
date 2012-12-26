//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>
#include <stdexcept>

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_container.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_statement.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>

#include <boost/range/iterator_range.hpp>

#include "cicada/neuron/convolution.hpp"

namespace cicada
{
  namespace neuron
  {
    Convolution::Convolution(size_type __frame, size_type __kW, size_type __dW)
      : frame(__frame), kW(__kW), dW(__dW)
    {
      weight = tensor_type(frame, 1);
      bias   = tensor_type(frame, 1);
    }
    
    Convolution::Convolution(const tensor_type& __weight, const tensor_type& __bias, size_type __kW, size_type __dW)
      : weight(__weight), bias(__bias), kW(__kW), dW(__dW)
    {
      frame = weight.rows();

      if (weight.rows() != bias.rows())
	throw std::runtime_error("invalid weight/bias");
      if (weight.cols() != bias.cols())
	throw std::runtime_error("invalid weight/bias");
      if (weight.cols() != 1)
	throw std::runtime_error("invalid weight");
      if (bias.cols() != 1)
	throw std::runtime_error("invalid bias");
    }

    void Convolution::forward(const tensor_type& data_input)
    {
      const size_type n_frame_input  = data_input.cols();
      const size_type n_frame_output = (n_frame_input - kW) / dW + 1;
      
      data_output.resize(frame, n_frame_output);
      
      for (size_type k = 0; k != n_frame_output; ++ k)
	data_output.col(k) = data_input.block(0, k * dW, data_input.rows(), kW).rowwise().sum().array() * weight.array() + bias.array();
    }
    
    void Convolution::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.resizeLike(data_input);
      gradient_input.setZero();
      
      const size_type n_frame_output = gradient_output.cols();
      
      for (size_type k = 0; k != n_frame_output; ++ k)
	gradient_input.block(0, k * dW, data_input.rows(), kW) += (gradient_output.col(k).transpose() * weight).replicate(1, kW);
    }

    void Convolution::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      const size_type n_frame_output = gradient_output.cols();

      gradient_weight.resizeLike(weight);
      gradient_bias.resizeLike(bias);

      gradient_weight.setZero();
      gradient_bias.setZero();
      
      for (size_type k = 0; k != n_frame_output; ++ k) {
	gradient_weight.array() += data_input.block(0, k * dW, data_input.rows(), kW).rowwise().sum().array() * gradient_output.col(k).array();
	gradient_bias   += gradient_output.col(k);
      }
    }

    template <typename Iterator>
    struct tensor_generator_grammar : boost::spirit::karma::grammar<Iterator, const Layer::tensor_type&()>
    {
      typedef Layer::tensor_type tensor_type;

      struct matrix_begin_func
      {
	template<class>
	struct result {
	  typedef const float* type;
	};
	
	const float* operator()(const Layer::tensor_type& tensor) const
	{
	  return tensor.data();
	}
      };

      struct rows_range_func
      {
	typedef boost::iterator_range<const float*> return_type;

	template<class, class>
	struct result {
	  typedef return_type type;
	};
	typedef const float* const_iterator;
	
	return_type operator()(const_iterator& value, int rows) const
	{
	  value += rows;
	  return return_type(value - rows, value);
	}
      };

      struct real_policy : boost::spirit::karma::real_policies<float>
      {
	static unsigned int precision(float)
	{
	  return 10;
	}
      };
      
      tensor_generator_grammar() : tensor_generator_grammar::base_type(tensor)
      {
	namespace karma = boost::spirit::karma;
	namespace standard = boost::spirit::standard;
        namespace phoenix = boost::phoenix;
	
	rows %= karma::lit('[') << (float10 % ',') << karma::lit(']');
	repeat_rows %= karma::lit(',') << rows;
	
	matrix %= (karma::lit('[')
		   << rows[karma::_1 = rows_range(karma::_r1, karma::_r2)]
		   << karma::repeat(karma::_r3 - 1)[karma::attr_cast(repeat_rows)[karma::_1 = rows_range(karma::_r1, karma::_r2)] ]
		   << karma::lit(']'));
	
	tensor = (karma::lit('{')
		  << karma::lit("\"row\"") << ':' << karma::int_[karma::_1 = phoenix::bind(&tensor_type::rows, karma::_val)]
	          << karma::lit(',')
		  << karma::lit("\"column\"") << ':' << karma::int_[karma::_1 = phoenix::bind(&tensor_type::cols, karma::_val)]
		  << karma::lit(',')
		  << karma::lit("\"tensor\"") << ':'
		  << matrix(matrix_begin(karma::_val),
			    phoenix::bind(&tensor_type::rows, karma::_val),
			    phoenix::bind(&tensor_type::cols, karma::_val))
		  << karma::lit('}'));
      }
      
      boost::spirit::karma::real_generator<float, real_policy> float10;
      
      boost::phoenix::function<rows_range_func> const rows_range;
      boost::phoenix::function<matrix_begin_func> const matrix_begin;
      
      boost::spirit::karma::rule<Iterator, boost::iterator_range<const float*>() > rows;
      boost::spirit::karma::rule<Iterator, boost::iterator_range<const float*>() > repeat_rows;
      boost::spirit::karma::rule<Iterator, void(const float*, int, int), boost::spirit::karma::locals<int /*_a*/> > matrix;
      boost::spirit::karma::rule<Iterator, const tensor_type&() > tensor;
    };

    std::ostream& Convolution::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;

      tensor_generator_grammar<std::ostream_iterator<char> > tensor;
      
      karma::generate(std::ostream_iterator<char>(os),
		      karma::lit('{')
		      << karma::lit("\"neuron\"") << ':' << karma::lit("\"convolution\"")
		      << ',' << karma::lit("\"weight\"") << ':' << tensor
		      << ',' << karma::lit("\"bias\"") << ':' << tensor
		      << ',' << karma::lit("\"kW\"") << ':' << karma::int_
		      << ',' << karma::lit("\"dW\"") << ':' << karma::int_
		      << '}',
		      weight,
		      bias,
		      kW,
		      dW);
      
      return os;
    }
  }
}
