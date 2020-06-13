//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>

#include <boost/spirit/include/karma.hpp>

#include "cicada/neuron/sequential.hpp"

#include <memory>

namespace cicada
{
  namespace neuron
  {
    
    void Sequential::forward(const tensor_type& data_input)
    {
      const tensor_type* data_curr = &data_input;
      
      layer_set_type::iterator liter_end = layers.end();
      for (layer_set_type::iterator liter = layers.begin(); liter != liter_end; ++ liter) {
	(*liter)->forward(*data_curr);
	data_curr = &(*liter)->data_output;
      }
      
      data_output = *data_curr;
    }
    
    void Sequential::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      const tensor_type* gradient_curr = &gradient_output;
      
      if (layers.size() > 1)
	for (size_type i = layers.size() - 1; i != 0; -- i) {
	  layers[i]->backward(layers[i - 1]->data_output, *gradient_curr);
	  
	  gradient_curr = &layers[i]->gradient_input;
	}
      
      layers.front()->backward(data_input, *gradient_curr);
      gradient_input = layers.front()->gradient_input;
    }

    void Sequential::accumulate(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      const tensor_type* gradient_curr = &gradient_output;
      
      if (layers.size() > 1)
	for (size_type i = layers.size() - 1; i != 0; -- i) {
	  layers[i]->accumulate(layers[i - 1]->data_output, *gradient_curr);
	  
	  gradient_curr = &layers[i]->gradient_input;
	}
      
      layers.front()->accumulate(data_input, *gradient_curr);
    }

    Sequential::layer_ptr_type Sequential::clone(const bool share) const
    {
      std::unique_ptr<Sequential> cloned(new Sequential(*this));
      
      for (size_type i = 0; i != layers.size(); ++ i)
	cloned->layers[i] = layers[i]->clone(share);
      
      return layer_ptr_type(cloned.release());
    }
    
    void Sequential::share(const layer_ptr_type& x)
    {
      if (! x)
	throw std::runtime_error("no layer?");
      
      const Sequential* other = dynamic_cast<const Sequential*>(x.get());
      
      if (! other || layers.size() != other->layers.size())
	throw std::runtime_error("invalid parameter sharing");
      
      for (size_type i = 0; i != layers.size(); ++ i)
	layers[i]->share(other->layers[i]);
    }
    
    std::ostream& Sequential::write(std::ostream& os) const
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      karma::generate(std::ostream_iterator<char>(os),
		      karma::lit('{')
		      << karma::lit("\"neuron\"") << ':' << karma::lit("\"sequential\"")
		      << ',' << karma::lit("\"layers\"") << ':' << (karma::lit('[') << (karma::stream % ',') << karma::lit(']'))
		      << '}',
		      layers);
      
      return os;
    }
  }
}
