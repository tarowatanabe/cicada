//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/apply.hpp>
#include <cicada/semiring.hpp>

#include <cicada/operation/apply.hpp>
#include <cicada/operation/functional.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    Apply::Apply(const std::string& parameter,
		 const model_type& __model,
		 const int __debug)
      : model(__model), weights(0), size(200), weights_one(false), exact(false), prune(false), grow(false), incremental(false), forced(false), debug(__debug)
    {
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "apply")
	throw std::runtime_error("this is not a feature-functin applier");

      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "size")
	  size = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "exact")
	  exact = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "prune")
	  prune = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "grow")
	  grow = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "incremental")
	  incremental = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "forced")
	  forced = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "weights")
	  weights = &base_type::weights(piter->second);
	else if (utils::ipiece(piter->first) == "weights-one")
	  weights_one = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "feature" || utils::ipiece(piter->first) == "feature-function")
	  model_local.push_back(feature_function_type::create(piter->second));
	else
	  std::cerr << "WARNING: unsupported parameter for apply: " << piter->first << "=" << piter->second << std::endl;
      }

      // default to prune...
      switch (int(exact) + prune + grow + incremental) {
      case 0: prune = true; break; // default to cube-prune
      case 1: break; // OK
      default:
	throw std::runtime_error("specify one of exact/prune/grow/incremental");
      }
    
      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
    }

    void Apply::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      typedef cicada::semiring::Logprob<double> weight_type;

      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type applied;

      model_type& __model = const_cast<model_type&>(! model_local.empty() ? model_local : model);
    
      // assignment...
      __model.assign(data.id, data.hypergraph, data.lattice, data.spans, data.targets, data.ngram_counts);
    
      if (forced)
	__model.apply_feature(true);
    
      weight_set_type weights_zero;
      const weight_set_type* weights_apply = (weights ? weights : &weights_zero);
    
      if (debug)
	std::cerr << "apply features: " << (exact ? "exact" : (incremental ? "incremental" : (grow ? "grow" : "prune"))) << std::endl;
    
      utils::resource start;
    
      // apply...
      if (exact)
	cicada::apply_exact(__model, hypergraph, applied);
      else if (incremental) {
	if (weights_one)
	  cicada::apply_incremental(__model, hypergraph, applied, weight_function_one<weight_type>(), size);
	else
	  cicada::apply_incremental(__model, hypergraph, applied, weight_function<weight_type>(*weights_apply), size);
      } else if (grow) {
	if (weights_one)
	  cicada::apply_cube_grow(__model, hypergraph, applied, weight_function_one<weight_type>(), size);
	else
	  cicada::apply_cube_grow(__model, hypergraph, applied, weight_function<weight_type>(*weights_apply), size);
      } else {
	if (weights_one)
	  cicada::apply_cube_prune(__model, hypergraph, applied, weight_function_one<weight_type>(), size);
	else
	  cicada::apply_cube_prune(__model, hypergraph, applied, weight_function<weight_type>(*weights_apply), size);
      }
    
      utils::resource end;
    
      __model.apply_feature(false);
    
      if (debug)
	std::cerr << "apply cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;

      if (debug)
	std::cerr << "# of nodes: " << applied.nodes.size()
		  << " # of edges: " << applied.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(applied.is_valid())
		  << std::endl;
	
      hypergraph.swap(applied);
    }

    void Apply::assign(const weight_set_type& __weights)
    {
      weights = &__weights;
    }
  };
};
