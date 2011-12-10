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
      : model(__model), weights(0), weights_assigned(0), size(200), weights_one(false), weights_fixed(false), exact(false), prune(false), grow(false), incremental(false), forced(false), sparse(false), dense(false), debug(__debug)
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
	else if (utils::ipiece(piter->first) == "sparse")
	  sparse = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "dense")
	  dense = utils::lexical_cast<bool>(piter->second);
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

      if (sparse && dense)
	throw std::runtime_error("either sparse|dense|all");
      
      if (sparse) {
	model_type model_sparse;

	model_type& __model = const_cast<model_type&>(! model_local.empty() ? model_local : model);

	for (model_type::const_iterator iter = __model.begin(); iter != __model.end(); ++ iter)
	  if ((*iter)->sparse_feature())
	    model_sparse.push_back(*iter);
	
	if (model_sparse.empty())
	  throw std::runtime_error("we have no sparse features");
	
	model_local.swap(model_sparse);
      } else if (dense) {
	model_type model_dense;

	model_type& __model = const_cast<model_type&>(! model_local.empty() ? model_local : model);

	for (model_type::const_iterator iter = __model.begin(); iter != __model.end(); ++ iter)
	  if (! (*iter)->sparse_feature())
	    model_dense.push_back(*iter);
	
	if (model_dense.empty())
	  throw std::runtime_error("we have no dense features");
	
	model_local.swap(model_dense);
      }

      if (const_cast<model_type&>(! model_local.empty() ? model_local : model).empty())
	throw std::runtime_error("no features to apply?");
      
      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
      
      if (weights || weights_one)
	weights_fixed = true;
      
      if (! weights)
	weights = &base_type::weights();
      
      name = (std::string("apply-")
	      + std::string(sparse ? "sparse-" : (dense ? "dense-" : ""))
	      + (exact ? "exact" : (incremental ? "incremental" : (grow ? "grow" : "prune"))));
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
      
      const weight_set_type* weights_apply = (weights_assigned ? weights_assigned : &(weights->weights));
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
    
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
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << " thread time: " << (end.thread_time() - start.thread_time())
		  << std::endl;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << applied.nodes.size()
		  << " # of edges: " << applied.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(applied.is_valid())
		  << std::endl;
      
      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += applied.nodes.size();
      stat.edge += applied.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
      
      hypergraph.swap(applied);
    }

    void Apply::assign(const weight_set_type& __weights)
    {
      if (! weights_fixed)
	weights_assigned = &__weights;
    }
  };
};
