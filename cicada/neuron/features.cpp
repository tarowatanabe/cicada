//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <iterator>
#include <stdexcept>
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

#include "cicada/neuron/features.hpp"
#include "cicada/feature_vector.hpp"

#include "utils/json_string_generator.hpp"

namespace cicada
{
  namespace neuron
  {
    void Features::forward(const tensor_type& data_input)
    {
      typedef cicada::FeatureVector<double> feature_vector_type;

      if (data_input.cols() != 1)
	throw std::runtime_error("invalid input");
      
      data_output.resize(features.size(), 1);
      data_output.setZero();
      
      const feature_vector_type& feats = *static_cast<const feature_vector_type*>((void*) &data_input);
      
      feature_vector_type::const_iterator fiter_end = feats.end();
      for (feature_vector_type::const_iterator fiter = feats.begin(); fiter != fiter_end; ++ fiter) {
	feature_map_type::const_iterator iter = features.find(fiter->first);
	
	if (iter != features.end())
	  data_output.col(0)[iter - features.begin()] = fiter->second;
      }
    }
    
    void Features::backward(const tensor_type& data_input, const tensor_type& gradient_output)
    {
      typedef cicada::FeatureVector<double> feature_vector_type;

      if (data_input.cols() != 1)
	throw std::runtime_error("invalid input");
      if (gradient_output.cols() != 1)
	throw std::runtime_error("invalid gradient output");
      
      const feature_vector_type& feats = *static_cast<const feature_vector_type*>((void*) &data_input);
      
      gradient_input.resize(features.size(), 1);
      gradient_input.setZero();
      
      feature_vector_type::const_iterator fiter_end = feats.end();
      for (feature_vector_type::const_iterator fiter = feats.begin(); fiter != fiter_end; ++ fiter) {
	feature_map_type::const_iterator iter = features.find(fiter->first);
	
	if (iter == features.end()) continue;
	
	if (gradient_input.rows() >= fiter->first.id()) {
	  const size_type feats_first = gradient_input.rows();
	  const size_type feats_last  = fiter->first.id() + 1;
	  
	  gradient_input.conservativeResize(fiter->first.id() + 1, 1);
	  gradient_input.block(feats_first, 0, feats_last - feats_first, 1).setZero();
	}
	
	gradient_input.col(0)[fiter->first.id()] = gradient_output.col(0)[iter - features.begin()];
      }
    }
    
    std::ostream& Features::write(std::ostream& os) const
    {
      typedef std::vector<feature_type, std::allocator<feature_type> > feats_type;

      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;

      utils::json_string_generator<std::ostream_iterator<char> > feature;
      
      karma::generate(std::ostream_iterator<char>(os),
		      karma::lit('{')
		      << karma::lit("\"neuron\"") << ':' << karma::lit("\"features\"")
		      << ',' << karma::lit("\"features\"") << ':' << karma::lit('[') << (feature % ',') << karma::lit(']')
		      << '}',
		      feats_type(features.begin(), features.end()));
      
      return os;
    }
  }
}
