//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/push_bos_eos.hpp>

#include <cicada/operation/push_bos_eos.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    PushBosEos::PushBosEos(const std::string& parameter, const int __debug)
      :  base_type("push-bos/eos"), debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "push-bos-eos")
	throw std::runtime_error("this is not a bos-eos pusher");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	std::cerr << "WARNING: unsupported parameter for remove-annotation: " << piter->first << "=" << piter->second << std::endl;
    }

    void PushBosEos::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type pushed;
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      utils::resource start;
	
      cicada::push_bos_eos(hypergraph, pushed);
	
      utils::resource end;
	
      if (debug)
	std::cerr << "push bos/eos cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
	
      if (debug)
	std::cerr << "push bos/eos: " << data.id
		  << " # of nodes: " << pushed.nodes.size()
		  << " # of edges: " << pushed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(pushed.is_valid())
		  << std::endl;
      
      hypergraph.swap(pushed);
    }
  };
};
