//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>

#include <cicada/operation/remove_head.hpp>
#include <cicada/remove_head.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    RemoveHead::RemoveHead(const std::string& parameter, const int __debug)
      :  base_type("remove-head"),
	 debug(__debug)
    {
      typedef cicada::Parameter param_type;
      
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "remove-head")
	throw std::runtime_error("this is not an head remover");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	std::cerr << "WARNING: unsupported parameter for remove-head: " << piter->first << "=" << piter->second << std::endl;
      }
    }

    void RemoveHead::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      hypergraph_type& hypergraph = data.hypergraph;
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      utils::resource start;
      
      hypergraph_type removed;
      cicada::remove_head(hypergraph, removed);
      
      utils::resource end;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << " thread time: " << (end.thread_time() - start.thread_time())
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
      stat.thread_time  += (end.thread_time() - start.thread_time());
      
      hypergraph.swap(removed);
    }
  };
};
