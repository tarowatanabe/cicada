//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>

#include <iostream>

#include <cicada/parameter.hpp>

#include <cicada/operation/functional.hpp>
#include <cicada/operation/posterior.hpp>
#include <cicada/posterior.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    Posterior::Posterior(const std::string& parameter, const int __debug)
      :  weights(0), weights_assigned(0), scale(1.0),
	 weights_one(false), weights_fixed(false), weights_extra(),
	 semiring_tropical(false), semiring_logprob(false), semiring_log(false),
	 debug(__debug)
    {
      typedef cicada::Parameter param_type;
      
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "posterior")
	throw std::runtime_error("this is not a posterior computer");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else if (utils::ipiece(piter->first) == "weights")
	  weights = &base_type::weights(piter->second);
	else if (utils::ipiece(piter->first) == "weights-one")
	  weights_one = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "scale")
	  scale = utils::lexical_cast<double>(piter->second);
	else if (utils::ipiece(piter->first) == "semiring") {
	  const utils::ipiece name = piter->second;
	  
	  if (name == "tropical")
	    semiring_tropical = true;
	  else if (name == "logprob")
	    semiring_logprob = true;
	  else if (name == "log")
	    semiring_log = true;
	  else
	    throw std::runtime_error("unknown semiring: " + piter->second);
	} else if (utils::ipiece(piter->first) == "weight") {
	  namespace qi = boost::spirit::qi;
	  namespace standard = boost::spirit::standard;

	  std::string::const_iterator iter = piter->second.begin();
	  std::string::const_iterator iter_end = piter->second.end();

	  std::string name;
	  double      value;
	  
	  if (! qi::phrase_parse(iter, iter_end,
				 qi::lexeme[+(!(qi::lit('=') >> qi::double_ >> (standard::space | qi::eoi))
					      >> (standard::char_ - standard::space))]
				 >> '='
				 >> qi::double_,
				 standard::blank, name, value) || iter != iter_end)
	    throw std::runtime_error("weight parameter parsing failed");
	  
	  weights_extra[name] = value;
	} else
	  std::cerr << "WARNING: unsupported parameter for posterior: " << piter->first << "=" << piter->second << std::endl;
      }

      if (int(semiring_tropical) + semiring_logprob + semiring_log > 1)
	throw std::runtime_error("you can specify one of tropical, logprob, log");
      
      if (int(semiring_tropical) + semiring_logprob + semiring_log == 0)
	semiring_tropical = true;
      
      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");

      if (weights_one && ! weights_extra.empty())
	throw std::runtime_error("you have extra weights, but specified all-one parameter");
      
      if (weights || weights_one)
	weights_fixed = true;
      
      if (! weights)
	weights = &base_type::weights();

      if (name.empty())
	name = "posterior";
      
      // assign feature-name!
      feature_name = static_cast<const std::string&>(name);
    }

    void Posterior::assign(const weight_set_type& __weights)
    {
      if (! weights_fixed)
	weights_assigned = &__weights;
    }

    void Posterior::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      hypergraph_type& hypergraph = data.hypergraph;
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;

      const weight_set_type* weights_posterior = (weights_assigned ? weights_assigned : &(weights->weights));

      hypergraph_type computed;
      
      utils::resource start;

      if (weights_one) {
	if (semiring_tropical)
	  cicada::posterior(hypergraph, computed, weight_scaled_function_one<cicada::semiring::Tropical<double> >(scale), feature_name);
	else if (semiring_logprob)
	  cicada::posterior(hypergraph, computed, weight_scaled_function_one<cicada::semiring::Logprob<double> >(scale), feature_name);
	else
	  cicada::posterior(hypergraph, computed, weight_scaled_function_one<cicada::semiring::Log<double> >(scale), feature_name);
      } else if (! weights_extra.empty()) {
	if (semiring_tropical)
	  cicada::posterior(hypergraph, computed, weight_scaled_function_extra<cicada::semiring::Tropical<double> >(*weights_posterior, scale, weights_extra.begin(), weights_extra.end()), feature_name);
	else if (semiring_logprob)
	  cicada::posterior(hypergraph, computed, weight_scaled_function_extra<cicada::semiring::Logprob<double> >(*weights_posterior, scale, weights_extra.begin(), weights_extra.end()), feature_name);
	else
	  cicada::posterior(hypergraph, computed, weight_scaled_function_extra<cicada::semiring::Log<double> >(*weights_posterior, scale, weights_extra.begin(), weights_extra.end()), feature_name);
      } else {
	if (semiring_tropical)
	  cicada::posterior(hypergraph, computed, weight_scaled_function<cicada::semiring::Tropical<double> >(*weights_posterior, scale), feature_name);
	else if (semiring_logprob)
	  cicada::posterior(hypergraph, computed, weight_scaled_function<cicada::semiring::Logprob<double> >(*weights_posterior, scale), feature_name);
	else
	  cicada::posterior(hypergraph, computed, weight_scaled_function<cicada::semiring::Log<double> >(*weights_posterior, scale), feature_name);
      }
      
      utils::resource end;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << " thread time: " << (end.thread_time() - start.thread_time())
		  << std::endl;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << computed.nodes.size()
		  << " # of edges: " << computed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(computed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += computed.nodes.size();
      stat.edge += computed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
      
      hypergraph.swap(computed);
    }
  };
};
