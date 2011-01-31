//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/expected_ngram.hpp>

#include <cicada/operation/expected_ngram.hpp>
#include <cicada/operation/functional.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    ExpectedNGram::ExpectedNGram(const std::string& parameter, const int __debug)
      : order(0), bos_eos(false), weights(0), weights_one(false), scale(1.0), debug(__debug)
    {
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "expected-ngram")
	throw std::runtime_error("this is not an expected-ngram computer"); 
    
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = boost::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "bos-eos")
	  bos_eos = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "weights")
	  weights = &base_type::weights(piter->second);
	else if (utils::ipiece(piter->first) == "weights-one")
	  weights_one = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "scale")
	  scale = boost::lexical_cast<double>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
      }

      if (order <= 0)
	throw std::runtime_error("order must be positive");
    
      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
    }

    void ExpectedNGram::operator()(data_type& data) const
    {
      typedef cicada::semiring::Logprob<double> weight_type;
      
      ngram_count_set_type& ngram_counts = data.ngram_counts;
      const hypergraph_type& hypergraph = data.hypergraph;
      
      ngram_counts.clear();
      if (! hypergraph.is_valid()) return;
      
      weight_set_type weights_zero;
      const weight_set_type* weights_apply = (weights ? weights : &weights_zero);
          
      utils::resource ngram_start;
    
      if (weights_one)
	cicada::expected_ngram(hypergraph, weight_scaled_function_one<weight_type>(scale), ngram_counts, order, bos_eos);
      else
	cicada::expected_ngram(hypergraph, weight_scaled_function<weight_type>(*weights_apply, scale), ngram_counts, order, bos_eos);

      utils::resource ngram_end;
    
      if (debug)
	std::cerr << "expected ngram cpu time: " << (ngram_end.cpu_time() - ngram_start.cpu_time())
		  << " user time: " << (ngram_end.user_time() - ngram_start.user_time())
		  << std::endl;
	
      if (debug)
	std::cerr << "expected ngram counts: size: " << ngram_counts.size() << std::endl;
    }

    void ExpectedNGram::assign(const weight_set_type& __weights)
    {
      weights = &__weights;
    }
    
  };
};
