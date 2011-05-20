//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "cicada/operation/functional.hpp"
#include "cicada/operation/traversal.hpp"
#include "cicada/operation/viterbi.hpp"

#include <cicada/parameter.hpp>
#include <cicada/viterbi.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>
#include <utils/hashmurmur.hpp>

#include <google/dense_hash_map>

namespace cicada
{
  namespace operation
  {
    Viterbi::Viterbi(const std::string& parameter, const int __debug)
      : weights(0), weights_assigned(0), weights_one(false), weights_fixed(false),
	semiring_tropical(false), semiring_logprob(false), semiring_log(false),
	debug(__debug)
    {
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "viterbi")
	throw std::runtime_error("this is not a viterbir");

      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "weights")
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
	  std::cerr << "WARNING: unsupported parameter for viterbi: " << piter->first << "=" << piter->second << std::endl;
      }
      
      
      if (int(semiring_tropical) + semiring_logprob + semiring_log == 0)
	semiring_tropical = true;
    
      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
      
      if (weights || weights_one)
	weights_fixed = true;
      
      if (! weights)
	weights = &base_type::weights();
      
      name = std::string("viterbi");
    }
    
    void Viterbi::operator()(data_type& data) const
    {
      typedef hypergraph_type::id_type id_type;
      typedef std::vector<id_type, std::allocator<id_type> > head_set_type;
      typedef google::dense_hash_map<id_type, id_type, utils::hashmurmur<size_t>, std::equal_to<id_type> > node_map_type;
      typedef cicada::operation::edge_traversal::edge_set_type edge_set_type;

      if (! data.hypergraph.is_valid()) return;
	
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type viterbi;
      
      const weight_set_type* weights_viterbi = (weights_assigned ? weights_assigned : &(weights->weights));
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
	
      utils::resource start;
      
      edge_set_type edges;
      
      if (weights_one) {
	if (semiring_tropical) {
	  cicada::semiring::Tropical<double> weight;
	  cicada::viterbi(hypergraph, edges, weight, cicada::operation::edge_traversal(), weight_function_one<cicada::semiring::Tropical<double> >());
	} else if (semiring_logprob) {
	  cicada::semiring::Logprob<double> weight;
	  cicada::viterbi(hypergraph, edges, weight, cicada::operation::edge_traversal(), weight_function_one<cicada::semiring::Logprob<double> >());
	} else {
	  cicada::semiring::Log<double> weight;
	  cicada::viterbi(hypergraph, edges, weight, cicada::operation::edge_traversal(), weight_function_one<cicada::semiring::Log<double> >());
	}
      } else {
	if (semiring_tropical) {
	  cicada::semiring::Tropical<double> weight;
	  cicada::viterbi(hypergraph, edges, weight, cicada::operation::edge_traversal(), weight_function<cicada::semiring::Tropical<double> >(*weights_viterbi));
	} else if (semiring_logprob) {
	  cicada::semiring::Logprob<double> weight;
	  cicada::viterbi(hypergraph, edges, weight, cicada::operation::edge_traversal(), weight_function<cicada::semiring::Logprob<double> >(*weights_viterbi));
	} else {
	  cicada::semiring::Log<double> weight;
	  cicada::viterbi(hypergraph, edges, weight, cicada::operation::edge_traversal(), weight_function<cicada::semiring::Log<double> >(*weights_viterbi));
	}
      }
      
      head_set_type heads;
      edge_set_type tails;
      node_map_type node_maps;
      node_maps.set_empty_key(id_type(-1));

      heads.reserve(edges.size());

      id_type node_id = 0;
      edge_set_type::const_iterator eiter_end = edges.end();
      for (edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter) {
	std::pair<node_map_type::iterator, bool> result = node_maps.insert(std::make_pair(hypergraph.edges[*eiter].head, node_id));
	
	heads.push_back(result.first->second);
	node_id += result.second;
      }
      
      for (id_type node = 0; node != node_id; ++ node)
	viterbi.add_node();
      
      id_type edge_id = 0;
      for (edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter, ++ edge_id) {
	const hypergraph_type::edge_type& edge = hypergraph.edges[*eiter];
	
	tails.clear();
	hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer) {
	  node_map_type::const_iterator niter = node_maps.find(*titer);
	  if (niter == node_maps.end())
	    throw std::runtime_error("no node?");
	  
	  tails.push_back(niter->second);
	}
	
	hypergraph_type::edge_type& edge_viterbi = viterbi.add_edge(tails.begin(), tails.end());
	
	edge_viterbi.rule       = edge.rule;
	edge_viterbi.features   = edge.features;
	edge_viterbi.attributes = edge.attributes;
	
	viterbi.connect_edge(edge_viterbi.id, heads[edge_id]);
      }
      
      node_map_type::const_iterator niter = node_maps.find(hypergraph.goal);
      if (niter == node_maps.end())
	throw std::runtime_error("did not reach goal?");
      
      viterbi.goal = niter->second;
      viterbi.topologically_sort();
      
      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << viterbi.nodes.size()
		  << " # of edges: " << viterbi.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(viterbi.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += viterbi.nodes.size();
      stat.edge += viterbi.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      
      hypergraph.swap(viterbi);
    }
    
    void Viterbi::assign(const weight_set_type& __weights)
    {
      if (! weights_fixed)
	weights_assigned = &__weights;
    }
  };
};
