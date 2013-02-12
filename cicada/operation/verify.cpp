//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>

#include <cicada/operation/verify.hpp>
#include <cicada/verify.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    Verify::Verify(const std::string& parameter, const int __debug)
      :  base_type("verify"),
	 debug(__debug)
    {
      typedef cicada::Parameter param_type;
      
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "verify")
	throw std::runtime_error("this is not a verifier");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	std::cerr << "WARNING: unsupported parameter for remove-unary: " << piter->first << "=" << piter->second << std::endl;
      }
    }

    void Verify::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      hypergraph_type& hypergraph = data.hypergraph;
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      utils::resource start;
      
      cicada::verify(hypergraph);
      
      utils::resource end;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << " thread time: " << (end.thread_time() - start.thread_time())
		  << std::endl;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << hypergraph.nodes.size()
		  << " # of edges: " << hypergraph.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(hypergraph.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += hypergraph.nodes.size();
      stat.edge += hypergraph.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
    }
  };
};
