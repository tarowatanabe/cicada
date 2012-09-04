//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/push_weights.hpp>

#include <cicada/operation/push_weights.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    PushWeights::PushWeights(const std::string& parameter, const int __debug)
      :  left(false), frontier(false), root(false), debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "push-weights")
	throw std::runtime_error("this is not a weight pusher");
      
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "direction") {
	  const utils::ipiece dir = piter->second;
	  
	  if (dir == "left")
	    left = true;
	  else if (dir == "frontier")
	    frontier = true;
	  else if (dir == "root")
	    root = true;
	  else
	    throw std::runtime_error("unuspported direction: " + parameter);
	} else
	  std::cerr << "WARNING: unsupported parameter for weight pusher: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (int(left) + frontier + root == 0)
	throw std::runtime_error("what direction for weight pushing? left, frontier or root");
      
      if (int(left) + frontier + root > 1)
	throw std::runtime_error("only single direction for weight pushing");
      
      name = std::string("push-weights-") + (left ? "left" : (frontier ? "frontier" : "root"));
    }

    void PushWeights::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type pushed;
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      utils::resource start;
      
      if (left)
	cicada::push_weights_left(hypergraph, pushed);
      else if (frontier)
	cicada::push_weights_frontier(hypergraph, pushed);
      else if (root)
	cicada::push_weights_root(hypergraph, pushed);
      else
	throw std::runtime_error("unsupported direction!");
	
      utils::resource end;
	
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << " thread time: " << (end.thread_time() - start.thread_time())
		  << std::endl;
	
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << pushed.nodes.size()
		  << " # of edges: " << pushed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(pushed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += pushed.nodes.size();
      stat.edge += pushed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
      
      hypergraph.swap(pushed);
    }
  };
};
