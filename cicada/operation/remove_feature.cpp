//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>

#include <cicada/operation/remove_feature.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    RemoveFeature::RemoveFeature(const std::string& parameter, const int __debug)
      :  base_type("remove-feature"),
	 debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "remove-feature")
	throw std::runtime_error("this is not an feature remover");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "feature")
	  removes.push_back(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for remove-feature: " << piter->first << "=" << piter->second << std::endl;
      }

      if (removes.empty())
	throw std::runtime_error("no feature to remove?");

      std::sort(removes.begin(), removes.end());
    }

    void RemoveFeature::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      hypergraph_type& hypergraph = data.hypergraph;
      
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      utils::resource start;
      
      hypergraph_type removed = hypergraph;
      
      hypergraph_type::edge_set_type::iterator eiter_end = removed.edges.end();
      for (hypergraph_type::edge_set_type::iterator eiter = removed.edges.begin(); eiter != eiter_end; ++ eiter) {
	hypergraph_type::edge_type& edge = *eiter;
	
	remove_set_type::const_iterator riter_end = removes.end();
	for (remove_set_type::const_iterator riter = removes.begin(); riter != riter_end; ++ riter)
	  edge.features.erase(*riter);
      }
      
      utils::resource end;
	
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
	
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << removed.nodes.size()
		  << " # of edges: " << removed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(removed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += removed.nodes.size();
      stat.edge += removed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      
      hypergraph.swap(removed);
    }
  };
};
