//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "cicada/operation/functional.hpp"
#include "cicada/operation/prune.hpp"

#include <cicada/parameter.hpp>
#include <cicada/prune.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    Prune::Prune(const std::string& parameter, const int __debug)
      : weights(0), weights_assigned(0), kbest(0), edge(0), beam(-1), density(0.0), scale(1.0), weights_one(false), weights_fixed(false),
	semiring_tropical(false), semiring_logprob(false), semiring_log(false),
	debug(__debug)
    {
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "prune")
	throw std::runtime_error("this is not a pruner");

      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "beam")
	  beam = utils::lexical_cast<double>(piter->second);
	else if (utils::ipiece(piter->first) == "kbest")
	  kbest = utils::lexical_cast<size_t>(piter->second);
	else if (utils::ipiece(piter->first) == "edge")
	  edge = utils::lexical_cast<size_t>(piter->second);
	else if (utils::ipiece(piter->first) == "density")
	  density = utils::lexical_cast<double>(piter->second);
	else if (utils::ipiece(piter->first) == "scale")
	  scale = utils::lexical_cast<double>(piter->second);
	else if (utils::ipiece(piter->first) == "weights")
	  weights = &base_type::weights(piter->second);
	else if (utils::ipiece(piter->first) == "weights-one")
	  weights_one = utils::lexical_cast<bool>(piter->second);
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
	
	} else
	  std::cerr << "WARNING: unsupported parameter for prune: " << piter->first << "=" << piter->second << std::endl;
      }

      const bool beam_mode = beam >= 0.0;
      const bool density_mode = density >= 1.0;
      const bool kbest_mode = kbest > 0;
      const bool edge_mode = edge > 0;
      
      if (int(beam_mode) + density_mode + kbest_mode + edge_mode > 1)
	throw std::runtime_error("you can specify one of kbest, beam and density pruning");
      
      if (int(beam_mode) + density_mode + kbest_mode + edge_mode == 0)
	throw std::runtime_error("you may want to specify either kbest, beam or density pruning");

      if (int(semiring_tropical) + semiring_logprob + semiring_log == 0)
	semiring_tropical = true;
    
      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
      
      if (weights || weights_one)
	weights_fixed = true;
      
      if (! weights)
	weights = &base_type::weights();
    }
    
    void Prune::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
	
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type pruned;
      
      const weight_set_type* weights_prune = (weights_assigned ? weights_assigned : &(weights->weights));
      
      const bool beam_mode = beam >= 0.0;
      const bool density_mode = density >= 1.0;
      const bool kbest_mode = kbest > 0;
      const bool edge_mode = edge > 0;
	
      if (debug)
	std::cerr << "prune " << (edge_mode ? "edge" : (kbest_mode ? "kbest" : (beam_mode ? "beam" : "density"))) << ": " << data.id << std::endl;
	
      utils::resource prune_start;
      
      if (weights_one) {
	if (edge_mode) {
	  if (semiring_tropical)
	    cicada::prune_edge(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Tropical<double> >(scale), edge);
	  else if (semiring_logprob)
	    cicada::prune_edge(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Logprob<double> >(scale), edge);
	  else
	    cicada::prune_edge(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Log<double> >(scale), edge);
	} else if (kbest_mode) {
	  if (semiring_tropical)
	    cicada::prune_kbest(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Tropical<double> >(scale), kbest);
	  else if (semiring_logprob)
	    cicada::prune_kbest(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Logprob<double> >(scale), kbest);
	  else
	    cicada::prune_kbest(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Log<double> >(scale), kbest);
	} else if (beam_mode) {
	  if (semiring_tropical)
	    cicada::prune_beam(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Tropical<double> >(scale), beam);
	  else if (semiring_logprob)
	    cicada::prune_beam(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Logprob<double> >(scale), beam);
	  else
	    cicada::prune_beam(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Log<double> >(scale), beam);
	} else if (density_mode) {
	  if (semiring_tropical)
	    cicada::prune_density(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Tropical<double> >(scale), density);
	  else if (semiring_logprob)
	    cicada::prune_density(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Logprob<double> >(scale), density);
	  else
	    cicada::prune_density(hypergraph, pruned, weight_scaled_function_one<cicada::semiring::Log<double> >(scale), density);
	} else
	  throw std::runtime_error("what pruning?");
      } else {
	if (edge_mode) {
	  if (semiring_tropical)
	    cicada::prune_edge(hypergraph, pruned, weight_scaled_function<cicada::semiring::Tropical<double> >(*weights_prune, scale), edge);
	  else if (semiring_logprob)
	    cicada::prune_edge(hypergraph, pruned, weight_scaled_function<cicada::semiring::Logprob<double> >(*weights_prune, scale), edge);
	  else
	    cicada::prune_edge(hypergraph, pruned, weight_scaled_function<cicada::semiring::Log<double> >(*weights_prune, scale), edge);
	} else if (kbest_mode) {
	  if (semiring_tropical)
	    cicada::prune_kbest(hypergraph, pruned, weight_scaled_function<cicada::semiring::Tropical<double> >(*weights_prune, scale), kbest);
	  else if (semiring_logprob)
	    cicada::prune_kbest(hypergraph, pruned, weight_scaled_function<cicada::semiring::Logprob<double> >(*weights_prune, scale), kbest);
	  else
	    cicada::prune_kbest(hypergraph, pruned, weight_scaled_function<cicada::semiring::Log<double> >(*weights_prune, scale), kbest);
	} else if (beam_mode) {
	  if (semiring_tropical)
	    cicada::prune_beam(hypergraph, pruned, weight_scaled_function<cicada::semiring::Tropical<double> >(*weights_prune, scale), beam);
	  else if (semiring_logprob)
	    cicada::prune_beam(hypergraph, pruned, weight_scaled_function<cicada::semiring::Logprob<double> >(*weights_prune, scale), beam);
	  else
	    cicada::prune_beam(hypergraph, pruned, weight_scaled_function<cicada::semiring::Log<double> >(*weights_prune, scale), beam);
	} else if (density_mode) {
	  if (semiring_tropical)
	    cicada::prune_density(hypergraph, pruned, weight_scaled_function<cicada::semiring::Tropical<double> >(*weights_prune, scale), density);
	  else if (semiring_logprob)
	    cicada::prune_density(hypergraph, pruned, weight_scaled_function<cicada::semiring::Logprob<double> >(*weights_prune, scale), density);
	  else
	    cicada::prune_density(hypergraph, pruned, weight_scaled_function<cicada::semiring::Log<double> >(*weights_prune, scale), density);
	} else
	  throw std::runtime_error("what pruning?");
      }
	
      utils::resource prune_end;
    
      if (debug)
	std::cerr << "prune cpu time: " << (prune_end.cpu_time() - prune_start.cpu_time())
		  << " user time: " << (prune_end.user_time() - prune_start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << "prune: " << data.id
		  << " # of nodes: " << pruned.nodes.size()
		  << " # of edges: " << pruned.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(pruned.is_valid())
		  << std::endl;
    
      hypergraph.swap(pruned);

    }
    
    void Prune::assign(const weight_set_type& __weights)
    {
      if (! weights_fixed)
	weights_assigned = &__weights;
    }
  };
};
