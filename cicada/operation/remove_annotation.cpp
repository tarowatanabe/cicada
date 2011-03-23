//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/remove_annotation.hpp>

#include <cicada/operation/remove_annotation.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    RemoveAnnotation::RemoveAnnotation(const std::string& parameter, const int __debug)
      :  debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "remove-annotation")
	throw std::runtime_error("this is not an annotation remover");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	std::cerr << "WARNING: unsupported parameter for remove-annotation: " << piter->first << "=" << piter->second << std::endl;
    }

    void RemoveAnnotation::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type removed;
      
      if (debug)
	std::cerr << "remove annotation: " << data.id << std::endl;
      
      utils::resource start;
	
      cicada::remove_annotation(hypergraph, removed);
	
      utils::resource end;
	
      if (debug)
	std::cerr << "remove annotation cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
	
      if (debug)
	std::cerr << "remove annotation: " << data.id
		  << " # of nodes: " << removed.nodes.size()
		  << " # of edges: " << removed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(removed.is_valid())
		  << std::endl;
      
      hypergraph.swap(removed);
    }
  };
};
