// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__INTERSECT__HPP__
#define __CICADA__OPERATION__INTERSECT__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/intersect.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

namespace cicada
{
  namespace operation
  {
    class Intersect : public Operation
    {
    public:
      Intersect(const std::string& parameter, const int __debug)
	: lattice_mode(false), target_mode(false), debug(__debug)
      {
	typedef cicada::Parameter param_type;
	
	param_type param(parameter);
	if (param.name() != "intersect")
	  throw std::runtime_error("this is not a intersector");

	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "lattice") == 0)
	    lattice_mode = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "target") == 0)
	    target_mode = utils::lexical_cast<bool>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for intersect: " << piter->first << "=" << piter->second << std::endl;
	}

	if (lattice_mode && target_mode)
	  throw std::runtime_error("either lattice or target");
	
	if (! lattice_mode && ! target_mode)
	  target_mode = true;
      }

      void operator()(data_type& data) const
      {
	const sentence_set_type& targets = data.targets;
	hypergraph_type& hypergraph = data.hypergraph;
	
	hypergraph_type intersected;
	
	utils::resource start;
	
	if (lattice_mode)
	  cicada::intersect(hypergraph, data.lattice, intersected);
	else {
	  if (targets.empty())
	    throw std::runtime_error("no target?");
	  
	  lattice_type target(targets.front());
	  
	  cicada::intersect(hypergraph, target, intersected);
	}
    
	utils::resource end;
    
	if (debug)
	  std::cerr << "intersect cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << std::endl;
	
	if (debug)
	  std::cerr << "# of nodes: " << intersected.nodes.size()
		    << " # of edges: " << intersected.edges.size()
		    << " valid? " << utils::lexical_cast<std::string>(intersected.is_valid())
		    << std::endl;
    
	hypergraph.swap(intersected);
      }

      bool lattice_mode;
      bool target_mode;
  
      int debug;
    };

  };
};


#endif
