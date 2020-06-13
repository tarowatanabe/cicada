//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/parallel.hpp"

#include <stdexcept>
#include <memory>

namespace cicada
{
  namespace neuron
  {
    void Parallel::forward(const tensor_type& data_input)
    {
      sizes.resize(layers.size());
      data_output.resize(0, 0);
      
      for (size_type i = 0; i != layers.size(); ++ i) {
	if (dimension_input)
	  layers[i]->forward(data_input.col(i));
	else
	  layers[i]->forward(data_input.row(i).transpose());
	
	sizes[i] = layers[i]->data_output.rows();
	
	if (dimension_output) {
	  if (data_output.cols()) {
	    if (data_output.rows() != static_cast<int>(sizes[i]))
	      throw std::runtime_error("invalid concat");
	    
	    data_output.conservativeResize(Eigen::NoChange, data_output.cols() + 1);
	    data_output.col(data_output.cols() - 1) = data_output.col(0);
	  } else
	    data_output = layers[i]->data_output.col(0);
	} else {
	  if (data_output.rows()) {
	    data_output.conservativeResize(data_output.rows() + sizes[i], 1);
	    data_output.block(data_output.rows() - sizes[i], 0, sizes[i], 1) = data_output.col(0);
	  } else
	    data_output = layers[i]->data_output.col(0);
	}
      }

    }
    
    void Parallel::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      gradient_input.resizeLike(data_input);
      gradient_input.setZero();
      
      if (dimension_output) {
	for (size_type i = 0; i != sizes.size(); ++ i) {
	  if (dimension_input) {
	    layers[i]->backward(data_input.col(i), gradient_output.col(i));
	    gradient_input.col(i) += layers[i]->gradient_input;
	  } else {
	    layers[i]->backward(data_input.row(i).transpose(), gradient_output.col(i));
	    gradient_input.row(i) += layers[i]->gradient_input.transpose();
	  }
	}
      } else {
	size_type offset = 0;
	for (size_type i = 0; i != sizes.size(); ++ i) {
	  if (dimension_input) {
	    layers[i]->backward(data_input.col(i), gradient_output.block(offset, 0, sizes[i], 1));
	    gradient_input.col(i) += layers[i]->gradient_input;
	  } else {
	    layers[i]->backward(data_input.row(i).transpose(), gradient_output.block(offset, 0, sizes[i], 1));
	    gradient_input.row(i) += layers[i]->gradient_input.transpose();
	  }
	  
	  offset += sizes[i];
	}
      }
    }

    void Parallel::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      if (dimension_output) {
	for (size_type i = 0; i != sizes.size(); ++ i) {
	  if (dimension_input)
	    layers[i]->accumulate(data_input.col(i), gradient_output.col(i));
	  else
	    layers[i]->accumulate(data_input.row(i).transpose(), gradient_output.col(i));
	}
      } else {
	size_type offset = 0;
	for (size_type i = 0; i != sizes.size(); ++ i) {
	  if (dimension_input)
	    layers[i]->accumulate(data_input.col(i), gradient_output.block(offset, 0, sizes[i], 1));
	  else
	    layers[i]->accumulate(data_input.row(i).transpose(), gradient_output.block(offset, 0, sizes[i], 1));
	  
	  offset += sizes[i];
	}
      }
    }

    Parallel::layer_ptr_type Parallel::clone(const bool share) const
    {
      std::unique_ptr<Parallel> cloned(new Parallel(*this));
      
      for (size_type i = 0; i != layers.size(); ++ i)
	cloned->layers[i] = layers[i]->clone(share);
      
      return layer_ptr_type(cloned.release());
    }

    void Parallel::share(const layer_ptr_type& x)
    {
      if (! x)
	throw std::runtime_error("no layer?");
      
      const Parallel* other = dynamic_cast<const Parallel*>(x.get());
      
      if (! other || layers.size() != other->layers.size())
	throw std::runtime_error("invalid parameter sharing");
      
      for (size_type i = 0; i != layers.size(); ++ i)
	layers[i]->share(other->layers[i]);
    }

    std::ostream& Parallel::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os),
		      karma::lit('{')
		      << karma::lit("\"neuron\"") << ':' << karma::lit("\"parallel\"")
		      << ',' << karma::lit("\"input\"") << ':' << karma::bool_
		      << ',' << karma::lit("\"output\"") << ':' << karma::bool_
		      << ',' << karma::lit("\"layers\"") << ':' << (karma::lit('[') << (karma::stream % ',') << karma::lit(']'))
		      << '}',
		      dimension_input,
		      dimension_output,
		      layers);
      
      return os;
    }
  }
}
