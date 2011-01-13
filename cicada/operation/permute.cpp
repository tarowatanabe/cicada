//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/permute.hpp>

#include <cicada/operation/permute.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

#include <boost/functional/hash.hpp>

namespace cicada
{
  namespace operation
  {
    Permute::Permute(const std::string& parameter, const int __debug)
      : excludes(), size(0), debug(__debug)
    {
      typedef cicada::Parameter param_type;

      excludes.set_empty_key(symbol_type());
    
      param_type param(parameter);
      if (param.name() != "permute")
	throw std::runtime_error("this is not a permuter");

      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "size") == 0)
	  size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "exclude") == 0)
	  excludes.insert(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for permute: " << piter->first << "=" << piter->second << std::endl;
      }
    }
    
    void Permute::operator()(data_type& data) const
    {
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type permuted;
    
      if (debug)
	std::cerr << "permute" << std::endl;
    
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
	std::cerr << "# of nodes: " << permuted.nodes.size()
		  << " # of edges: " << permuted.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(permuted.is_valid())
		  << std::endl;
    
      hypergraph.swap(permuted);
    }
    
  };
};
