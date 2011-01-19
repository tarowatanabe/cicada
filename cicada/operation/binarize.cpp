//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/binarize.hpp>

#include <cicada/operation/binarize.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

namespace cicada
{
  namespace operation
  {

    Binarize::Binarize(const std::string& parameter, const int __debug)
      : order(-1), left(false), right(false), debug(__debug)
    {
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (param.name() != "binarize")
	throw std::runtime_error("this is not a binarizer");

      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "order") == 0)
	  order = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "direction") == 0) {
	  const std::string& dir = piter->second;
	
	  if (strcasecmp(dir.c_str(), "left") == 0)
	    left = true;
	  else if (strcasecmp(dir.c_str(), "right") == 0)
	    right = true;
	  else
	    throw std::runtime_error("unuspported direction: " + parameter);
	} else
	  std::cerr << "WARNING: unsupported parameter for binarize: " << piter->first << "=" << piter->second << std::endl;
      }
    
      if (! left && ! right)
	throw std::runtime_error("what direction?");
      if (left && right)
	throw std::runtime_error("we do not binarization in both directions!");
    }

    void Binarize::operator()(data_type& data) const
    {
      hypergraph_type binarized;
    
      if (debug)
	std::cerr << "binarization" << std::endl;
    
      utils::resource start;
    
      if (left)
	cicada::binarize_left(data.hypergraph, binarized, order);
      else if (right)
	cicada::binarize_right(data.hypergraph, binarized, order);
    
      utils::resource end;
    
      if (debug)
	std::cerr << "binarize cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << "# of nodes: " << binarized.nodes.size()
		  << " # of edges: " << binarized.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(binarized.is_valid())
		  << std::endl;
	
      data.hypergraph.swap(binarized);
    }
    
  };
};
