// -*- mode: c++ -*-

#ifndef __CICADA__OPERATION__EXPECTED_NGRAM__HPP__
#define __CICADA__OPERATION__EXPECTED_NGRAM__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/expected_ngram.hpp>

#include <cicada/operation/functional.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

namespace cicada
{
  namespace operation
  {
    class ExpectedNGram : public Operation
    {
    public:
      ExpectedNGram(const std::string& parameter, const int __debug)
	: order(0), bos_eos(false), weights(0), weights_one(false), yield_source(false), yield_target(false)
      {
	typedef cicada::Parameter param_type;
    
	param_type param(parameter);
	if (param.name() != "expected-ngram")
	  throw std::runtime_error("this is not an expected-ngram computer"); 
    
	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "order") == 0)
	    order = boost::lexical_cast<int>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "bos-eos") == 0)
	    bos_eos = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "weights") == 0)
	    weights = &base_type::weights(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "weights-one") == 0)
	    weights_one = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	    const std::string& value = piter->second;
	
	    if (strcasecmp(value.c_str(), "source") == 0)
	      yield_source = true;
	    else if (strcasecmp(value.c_str(), "target") == 0)
	      yield_target = true;
	  } else
	    std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
	}

	if (order <= 0)
	  throw std::runtime_error("order must be positive");
    
	if (weights && weights_one)
	  throw std::runtime_error("you have weights, but specified all-one parameter");
    
	if (yield_source && yield_target)
	  throw std::runtime_error("do not specify both source and target yield");
    
	if (! yield_source && ! yield_target)
	  yield_target = true;
      }

      void operator()(data_type& data) const
      {
	typedef cicada::semiring::Logprob<double> weight_type;

	ngram_count_set_type& ngram_counts = data.ngram_counts;
	const hypergraph_type& hypergraph = data.hypergraph;

	weight_set_type weights_zero;
	const weight_set_type* weights_apply = (weights ? weights : &weights_zero);
    
	ngram_counts.clear();
	
	utils::resource ngram_start;
    
	if (weights_one)
	  cicada::expected_ngram(hypergraph, weight_function_one<weight_type>(), ngram_counts, order, bos_eos, yield_source);
	else
	  cicada::expected_ngram(hypergraph, weight_function<weight_type>(*weights_apply), ngram_counts, order, bos_eos, yield_source);

	utils::resource ngram_end;
    
	if (debug)
	  std::cerr << "expected ngram cpu time: " << (ngram_end.cpu_time() - ngram_start.cpu_time())
		    << " user time: " << (ngram_end.user_time() - ngram_start.user_time())
		    << std::endl;
	
	if (debug)
	  std::cerr << "expected ngram counts: size: " << ngram_counts.size() << std::endl;
      }
  
      int order;
      bool bos_eos;
  
      const weight_set_type* weights;
      bool weights_one;
      bool yield_source;
      bool yield_target;
  
      int debug;
    };

  };
};


#endif
