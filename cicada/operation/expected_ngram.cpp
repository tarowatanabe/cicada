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
      : base_type("expected-ngram"),
	order(0), bos_eos(false), weights(0), weights_assigned(0), weights_one(false),weights_fixed(false), scale(1.0), debug(__debug)
    {
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "expected-ngram")
	throw std::runtime_error("this is not an expected-ngram computer"); 
    
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "bos-eos")
	  bos_eos = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "weights")
	  weights = &base_type::weights(piter->second);
	else if (utils::ipiece(piter->first) == "weights-one")
	  weights_one = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "scale")
	  scale = utils::lexical_cast<double>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
      }

      if (order <= 0)
	throw std::runtime_error("order must be positive");
    
      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
      
      if (weights || weights_one)
	weights_fixed = true;

      if (! weights)
	weights = &base_type::weights();
    }

    void ExpectedNGram::operator()(data_type& data) const
    {
      typedef cicada::semiring::Logprob<double> weight_type;
      
      ngram_count_set_type& ngram_counts = data.ngram_counts;
      const hypergraph_type& hypergraph = data.hypergraph;
      
      ngram_counts.clear();
      if (! hypergraph.is_valid()) return;
      
      const weight_set_type* weights_apply = (weights_assigned ? weights_assigned : &(weights->weights));
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
          
      utils::resource start;
    
      if (weights_one)
	cicada::expected_ngram(hypergraph, weight_scaled_function_one<weight_type>(scale), ngram_counts, order, bos_eos);
      else
	cicada::expected_ngram(hypergraph, weight_scaled_function<weight_type>(*weights_apply, scale), ngram_counts, order, bos_eos);

      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << " thread time: " << (end.thread_time() - start.thread_time())
		  << std::endl;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " counts: size: " << ngram_counts.size() << std::endl;

      
      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += ngram_counts.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
    }

    void ExpectedNGram::assign(const weight_set_type& __weights)
    {
      if (! weights_fixed)
	weights_assigned = &__weights;
    }
    
  };
};
