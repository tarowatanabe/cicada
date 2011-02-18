//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/sort_tail.hpp>

#include <cicada/operation/sort_tail.hpp>
#include <cicada/operation/functional.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    SortTail::SortTail(const std::string& parameter, const int __debug)
      : debug(__debug)
    {
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "sort-tail")
	throw std::runtime_error("this is not a tail sorter");

      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	std::cerr << "WARNING: unsupported parameter for tail sorter: " << piter->first << "=" << piter->second << std::endl;
    }
    
    void SortTail::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type sorted;

      if (debug)
	std::cerr << "sort tail: " << data.id
		  << " # of nodes: " << hypergraph.nodes.size()
		  << " # of edges: " << hypergraph.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(hypergraph.is_valid())
		  << std::endl;
    
      utils::resource start;

      cicada::sort_tail(hypergraph, sorted);
	
      utils::resource end;
    
      if (debug)
	std::cerr << "sort tail cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << "sorted tail: " << data.id
		  << " # of nodes: " << sorted.nodes.size()
		  << " # of edges: " << sorted.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(sorted.is_valid())
		  << std::endl;
    
      hypergraph.swap(sorted);
    }
  };
};
