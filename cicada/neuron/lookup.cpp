//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>
#include <stdexcept>
#include <memory>

#include <cstddef>
#include <stdint.h>

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_container.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_statement.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>

#include <boost/range/iterator_range.hpp>

#include "cicada/neuron/lookup.hpp"

namespace cicada
{
  namespace neuron
  {
    Lookup::Lookup(size_type __size) : size(__size), weight(new tensor_type()) {}
    
    Lookup::Lookup(const tensor_type& __weight) : size(__weight.rows()), weight(new tensor_type(__weight)) {}
    
    Lookup::Lookup(const tensor_ptr_type& __weight) : size(__weight->rows()), weight(__weight) {}

    void Lookup::forward(const tensor_type& data_input)
    {
      if (data_input.cols() != 1)
	throw std::runtime_error("invalid input");

      data_output.resize(size, data_input.rows());
      
      // we will cast data_input as a sequence of index...
      const uint32_t* first = reinterpret_cast<const uint32_t*>(data_input.data());
      const uint32_t* last = first + data_input.rows();
      for (size_type frame = 0; first != last; ++ first, ++ frame) {
	// perform resizing of weight
	if (*first >= weight->cols()) {
	  // from weight.cols() to *first + 1, assign random value...!
	  const size_type frame_first = weight->cols();
	  const size_type frame_last  = *first + 1;
	  
	  weight->conservativeResize(size, *first + 1);
	  for (size_type i = frame_first; i != frame_last; ++ i)
	    weight->col(i).setRandom();
	}
	
	data_output.col(frame) = weight->col(*first);
      }
    }
    
    void Lookup::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      // we will do nothing, since this must be the starting point!
      
    }
    
    void Lookup::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      if (data_input.cols() != 1)
	throw std::runtime_error("invalid input");
      
      if (gradient_output.rows() != size)
	throw std::runtime_error("invalid gradient output");
      
      gradient_weight.resizeLike(*weight);
      
      // we do not clear all, but only required part...
      //gradient_weight.setZero();
      {
	const uint32_t* first = reinterpret_cast<const uint32_t*>(data_input.data());
	const uint32_t* last = first + data_input.rows();
	for (/**/; first != last; ++ first)
	  gradient_weight.col(*first).setZero();
      }
      
      const uint32_t* first = reinterpret_cast<const uint32_t*>(data_input.data());
      const uint32_t* last = first + data_input.rows();
      for (size_type frame = 0; first != last; ++ first, ++ frame)
	gradient_weight.col(*first) += gradient_output.col(frame);
    }

    Lookup::layer_ptr_type Lookup::clone(const bool share) const
    {
      std::unique_ptr<Lookup> cloned(new Lookup(*this));
      
      if (! share)
	cloned->weight.reset(new tensor_type(*weight));
      
      return layer_ptr_type(cloned.release());
    }

    void Lookup::share(const layer_ptr_type& x)
    {
      if (! x)
	throw std::runtime_error("no layer?");
      
      const Lookup* other = dynamic_cast<const Lookup*>(x.get());
      
      if (! other)
	throw std::runtime_error("invalid parameter sharing");
      
      weight = other->weight;
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

    std::ostream& Lookup::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      tensor_generator_grammar<std::ostream_iterator<char> > tensor;

      karma::generate(std::ostream_iterator<char>(os),
		      karma::lit('{')
		      << karma::lit("\"neuron\"") << ':' << karma::lit("\"lookup\"")
		      << ',' << karma::lit("\"weight\"") << ':' << tensor
		      << '}',
		      *weight);
      
      return os;
    }
  }
}
