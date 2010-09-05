// -*- mode: c++ -*-

#ifndef __CICADA__OPERATION__GENERATE__HPP__
#define __CICADA__OPERATION__GENERATE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/generate.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

namespace cicada
{
  namespace operation
  {
    class GenerateEarley : public Operation
    {
    public:
      GenerateEarley(const grammar_type& __grammar,
		     const std::string& __goal,
		     const std::string& __non_terminal,
		     const bool __insertion,
		     const bool __deletion,
		     const int __debug)
	: debug(__debug)
      { }
  
      void operator()(data_type& data) const
      {
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type generated;
    
	if (debug)
	  std::cerr << "generation: earley" << std::endl;
    
	utils::resource start;
    
	cicada::generate_earley(hypergraph, generated);
    
	utils::resource end;
    
	if (debug)
	  std::cerr << "generate cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << std::endl;
    
	if (debug)
	  std::cerr << "generated: earley"
		    << " previous: # of nodes: " << hypergraph.nodes.size()
		    << " # of edges: " << hypergraph.edges.size()
		    << " generated: # of nodes: " << generated.nodes.size()
		    << " # of edges: " << generated.edges.size()
		    << " valid? " << utils::lexical_cast<std::string>(generated.is_valid())
		    << std::endl;
    
	hypergraph.swap(generated);
      }
  
      int debug;
    };

  };
};


#endif
