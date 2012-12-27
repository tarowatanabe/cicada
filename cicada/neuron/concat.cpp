//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>
#include <stdexcept>
#include <memory>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/concat.hpp"


namespace cicada
{
  namespace neuron
  {
    void Concat::forward(const tensor_type& data_input)
    {
      sizes.resize(layers.size());
      data_output.resize(0, 0);

      if (dimension) {
	size_set_type::iterator siter = sizes.begin();
	layer_set_type::iterator liter_end = layers.end();
	for (layer_set_type::iterator liter = layers.begin(); liter != liter_end; ++ liter, ++ siter) {
	  (*liter)->forward(data_input);
	  
	  *siter = (*liter)->data_output.rows();
	  
	  if (data_output.cols()) {
	    if (data_output.rows() != static_cast<int>(*siter))
	      throw std::runtime_error("invalid concat");

	    data_output.conservativeResize(Eigen::NoChange, data_output.cols() + 1);
	    data_output.col(data_output.cols() - 1) = data_output.col(0);
	  } else
	    data_output = (*liter)->data_output.col(0);
	}
      } else {
	size_set_type::iterator siter = sizes.begin();
	layer_set_type::iterator liter_end = layers.end();
	for (layer_set_type::iterator liter = layers.begin(); liter != liter_end; ++ liter, ++ siter) {
	  (*liter)->forward(data_input);
	  
	  *siter = (*liter)->data_output.rows();
	  
	  if (data_output.rows()) {
	    data_output.conservativeResize(data_output.rows() + *siter, 1);
	    data_output.block(data_output.rows() - *siter, 0, *siter, 1) = data_output.col(0);
	  } else
	    data_output = (*liter)->data_output.col(0);
	}
      }
    }
    
    void Concat::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.resizeLike(data_input);
      gradient_input.setZero();
      
      if (dimension) {
	for (size_type i = 0; i != sizes.size(); ++ i) {
	  layers[i]->backward(data_input, gradient_output.col(i));
	  gradient_input += layers[i]->gradient_input;
	}
      } else {
	size_type offset = 0;
	for (size_type i = 0; i != sizes.size(); ++ i) {
	  layers[i]->backward(data_input, gradient_output.block(offset, 0, sizes[i], 1));
	  gradient_input += layers[i]->gradient_input;
	  offset += sizes[i];
	}
      }
    }

    void Concat::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      if (dimension) {
	for (size_type i = 0; i != sizes.size(); ++ i)
	  layers[i]->accumulate(data_input, gradient_output.col(i));
      } else {
	size_type offset = 0;
	for (size_type i = 0; i != sizes.size(); ++ i) {
	  layers[i]->accumulate(data_input, gradient_output.block(offset, 0, sizes[i], 1));
	  offset += sizes[i];
	}
      }
    }

    Concat::layer_ptr_type Concat::clone(const bool share) const
    {
      std::auto_ptr<Concat> cloned(new Concat(*this));
      
      for (size_type i = 0; i != layers.size(); ++ i)
	cloned->layers[i] = layers[i]->clone(share);
      
      return layer_ptr_type(cloned.release());
    }

    void Concat::share(const layer_ptr_type& x)
    {
      if (! x)
	throw std::runtime_error("no layer?");
      
      const Concat* other = dynamic_cast<const Concat*>(x.get());
      
      if (! other || layers.size() != other->layers.size())
	throw std::runtime_error("invalid parameter sharing");
      
      for (size_type i = 0; i != layers.size(); ++ i)
	layers[i]->share(other->layers[i]);
    }

    std::ostream& Concat::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os),
		      karma::lit('{')
		      << karma::lit("\"neuron\"") << ':' << karma::lit("\"concat\"")
		      << ',' << karma::lit("\"dimension\"") << ':' << karma::bool_
		      << ',' << karma::lit("\"layers\"") << ':' << (karma::lit('[') << (karma::stream % ',') << karma::lit(']'))
		      << '}',
		      dimension,
		      layers);
      
      return os;
    }
  }
}
