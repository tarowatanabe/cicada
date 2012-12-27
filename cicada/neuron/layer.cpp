//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_container.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_statement.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>

#include <boost/thread.hpp>

#include <vector>

#include "cicada/neuron/neuron.hpp"

#include "utils/thread_specific_ptr.hpp"

namespace cicada
{
  namespace neuron
  {
    template <typename Iterator>
    struct layer_parser_grammar : boost::spirit::qi::grammar<Iterator, Layer::layer_ptr_type(), boost::spirit::standard::space_type>
    {
      typedef Layer::tensor_type tensor_type;

      typedef Layer::layer_ptr_type ptr_type;
      typedef std::vector<ptr_type, std::allocator<ptr_type> > layer_ptr_set_type;

      struct matrix_begin_func
      {
	template<class>
	struct result {
	  typedef float* type;
	};
	
	float* operator()(Layer::tensor_type& tensor) const
	{
	  return tensor.data();
	}
      };
      
      layer_parser_grammar() : layer_parser_grammar::base_type(layer)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	namespace phoenix = boost::phoenix;

	tensor = (qi::lit('{')
		  >> (qi::lit("\"row\"") >> ':' >> qi::int_ >> ','>> qi::lit("\"column\"") >> ':' >> qi::int_)
		  [qi::_val = phoenix::construct<Layer::tensor_type>(qi::_1, qi::_2)]
		  >> ',' >> qi::lit("\"tensor\"") >> ':' 
		  >> (qi::lit('[') [qi::_a = matrix_begin(qi::_val),
				    qi::_b = qi::_a + (phoenix::bind(&tensor_type::rows, qi::_val)
						       * phoenix::bind(&tensor_type::cols, qi::_val))]
		      >> ((qi::lit('[')
			   >> ((qi::float_ [qi::_pass = qi::_a != qi::_b, *qi::_a = qi::_1, ++ qi::_a]) % ',')
			   >> qi::lit(']')) % ',')
		      >> qi::lit(']'))
		  >> qi::lit('}'));
	
	layer = (qi::lit('{') >> qi::lit("\"neuron\"") >> ':'
		 >> (qi::lit("\"abs\"")
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Abs>())]
		     | (qi::lit("\"concat\"")
			>> ',' >> qi::lit("\"dimension\"") >> ':' >> qi::bool_
			>> ',' >> qi::lit("\"layers\"") >> ':' >> layers)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Concat>(phoenix::begin(qi::_2), phoenix::end(qi::_2), qi::_1))]
		     | (qi::lit("\"concat\"")
			>> ',' >> qi::lit("\"dimension\"") >> ':' >> qi::bool_)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Concat>(qi::_1))]
		     | (qi::lit("\"convolution\"")
			>> ',' >> qi::lit("\"frame\"") >> ':' >> qi::int_
			>> ',' >> qi::lit("\"kW\"") >> ':' >> qi::int_
			>> ',' >> qi::lit("\"dW\"") >> ':' >> qi::int_
			)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Convolution>(qi::_1, qi::_2, qi::_3))]
		     | (qi::lit("\"convolution\"")
			>> ',' >> qi::lit("\"weight\"") >> ':' >> tensor
			>> ',' >> qi::lit("\"bias\"") >> ':' >> tensor
			>> ',' >> qi::lit("\"kW\"") >> ':' >> qi::int_
			>> ',' >> qi::lit("\"dW\"") >> ':' >> qi::int_
			)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Convolution>(qi::_1, qi::_2, qi::_3, qi::_4))]
		     | qi::lit("\"copy\"")
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Copy>())]
		     | qi::lit("\"exp\"")
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Exp>())]
		     | (qi::lit("\"hardtanh\"") | qi::lit("\"hard-tanh\""))
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::HardTanh>())]
		     | (qi::lit("\"linear\"")
			>> ',' >> qi::lit("\"input\"") >> ':' >> qi::int_
			>> ',' >> qi::lit("\"output\"") >> ':' >> qi::int_)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Linear>(qi::_1, qi::_2))]
		     | (qi::lit("\"linear\"")
			>> ',' >> qi::lit("\"weight\"") >> ':' >> tensor
			>> ',' >> qi::lit("\"bias\"") >> ':' >> tensor)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Linear>(qi::_1, qi::_2))]
		     | qi::lit("\"log\"")
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Log>())]
		     | (qi::lit("\"lookup\"")
			>> ',' >> qi::lit("\"size\"") >> ':' >> qi::int_)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Lookup>(qi::_1))]
		     | (qi::lit("\"lookup\"")
			>> ',' >> qi::lit("\"weight\"") >> ':' >> tensor)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Lookup>(qi::_1))]
		     | (qi::lit("\"max\"")
			>> ',' >> qi::lit("\"dimension\"") >> ':' >> qi::bool_)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Max>(qi::_1))]
		     | (qi::lit("\"mean\"")
			>> ',' >> qi::lit("\"dimension\"") >> ':' >> qi::bool_)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Mean>(qi::_1))]
		     | (qi::lit("\"min\"")
			>> ',' >> qi::lit("\"dimension\"") >> ':' >> qi::bool_)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Min>(qi::_1))]
		     | (qi::lit("\"parallel\"")
			>> ',' >> qi::lit("\"input\"") >> ':' >> qi::bool_
			>> ',' >> qi::lit("\"output\"") >> ':' >> qi::bool_
			>> ',' >> qi::lit("\"layers\"") >> ':' >> layers)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Parallel>(phoenix::begin(qi::_3), phoenix::end(qi::_3), qi::_1, qi::_2))]
		     | (qi::lit("\"parallel\"")
			>> ',' >> qi::lit("\"input\"") >> ':' >> qi::bool_
			>> ',' >> qi::lit("\"output\"") >> ':' >> qi::bool_)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Parallel>(qi::_1, qi::_2))]
		     | (qi::lit("\"power\"")
			>> ',' >> qi::lit("\"power\"") >> ':' >> qi::double_)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Power>(qi::_1))]
		     | (qi::lit("\"sequential\"")
			>> ',' >> qi::lit("\"layers\"") >> ':' >> layers)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Sequential>(phoenix::begin(qi::_1), phoenix::end(qi::_1)))]
		     | qi::lit("\"sequential\"")
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Sequential>())]
		     | qi::lit("\"sigmoid\"")
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Sigmoid>())]
		     | qi::lit("\"softmax\"")
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::SoftMax>())]
		     | qi::lit("\"sqrt\"")
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Sqrt>())]
		     | qi::lit("\"square\"")
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Square>())]
		     | (qi::lit("\"sum\"")
			>> ',' >> qi::lit("\"dimension\"") >> ':' >> qi::bool_)
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Sum>(qi::_1))]
		     | qi::lit("\"tanh\"")
		     [qi::_val = phoenix::construct<ptr_type>(phoenix::new_<neuron::Tanh>())]
		     )
		 >> '}');
	
	layers %= qi::lit('[') >> (layer % ',') >> qi::lit(']');
      }
      
      typedef boost::spirit::standard::space_type space_type;

      boost::phoenix::function<matrix_begin_func> const matrix_begin;
      
      boost::spirit::qi::rule<Iterator, Layer::tensor_type(), space_type, boost::spirit::qi::locals<float* /*_a*/, float* /*_b*/> > tensor;
      boost::spirit::qi::rule<Iterator, layer_ptr_set_type(), space_type> layers;
      boost::spirit::qi::rule<Iterator, ptr_type(), space_type> layer;
    };
    
    Layer::layer_ptr_type Layer::construct(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;

      typedef layer_parser_grammar<std::string::const_iterator> parser_type;

      layer_ptr_type data;
      parser_type parser;
      if (qi::phrase_parse(iter, end, parser, standard::space, data))
	return data;
      else
	return layer_ptr_type();
    }
  
    Layer::layer_ptr_type Layer::construct(const utils::piece& data)
    {
      std::string::const_iterator iter(data.begin());
      std::string::const_iterator iter_end(data.end());
      
      layer_ptr_type layer = construct(iter, iter_end);
      
      if (iter != iter_end)
	throw std::runtime_error("neuron layer parsing failed");
      
      return layer;
    }
  }
}
