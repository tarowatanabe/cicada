//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/permute.hpp>

#include <cicada/operation/permute.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

#include <boost/functional/hash.hpp>

namespace cicada
{
  namespace operation
  {
    Permute::Permute(const std::string& parameter, const int __debug)
      : base_type("permute"),
	excludes(), size(0), debug(__debug)
    {
      typedef cicada::Parameter param_type;

      excludes.set_empty_key(symbol_type());
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "permute")
	throw std::runtime_error("this is not a permuter");

      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "size")
	  size = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "exclude")
	  excludes.insert(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for permute: " << piter->first << "=" << piter->second << std::endl;
      }
    }
    
    void Permute::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;

      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type permuted;
    
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
    
      utils::resource start;
	
      if (excludes.empty())
	cicada::permute(hypergraph, permuted, size);
      else
	cicada::permute(hypergraph, permuted, Filter(excludes), size);
	
      utils::resource end;
    
      if (debug)
	std::cerr << "permute cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << "permute: " << data.id 
		  << " # of nodes: " << permuted.nodes.size()
		  << " # of edges: " << permuted.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(permuted.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += permuted.nodes.size();
      stat.edge += permuted.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
    
      hypergraph.swap(permuted);
    }
    
  };
};
