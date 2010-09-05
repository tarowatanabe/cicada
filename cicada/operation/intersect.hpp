// -*- mode: c++ -*-

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
      Intersect(const int __debug)
	: debug(__debug) {}

      void operator()(data_type& data) const
      {
	const sentence_set_type& targets = data.targets;
	hypergraph_type& hypergraph = data.hypergraph;
    
	if (targets.empty())
	  throw std::runtime_error("no target?");
    
	lattice_type target(targets.front());

	hypergraph_type intersected;
    
	utils::resource start;
    
	cicada::intersect(hypergraph, target, intersected);
    
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
  
      int debug;
    };

  };
};


#endif
