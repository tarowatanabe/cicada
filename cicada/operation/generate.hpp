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
      GenerateEarley(const std::string& parameter,
		     const grammar_type& __grammar,
		     const std::string& __goal,
		     const std::string& __non_terminal,
		     const bool __insertion,
		     const bool __deletion,
		     const int __debug)
	: depth(0), width(0), context(false), debug(__debug)
      { 
	typedef cicada::Parameter param_type;
    
	param_type param(parameter);
	if (param.name() != "generate-earley")
	  throw std::runtime_error("this is not a Earley generator");
	
	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "depth") == 0)
	    depth = boost::lexical_cast<int>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "width") == 0)
	    width = boost::lexical_cast<int>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "context") == 0)
	    context = utils::lexical_cast<bool>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for generator: " << piter->first << "=" << piter->second << std::endl;
	}
      }
  
      void operator()(data_type& data) const
      {
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type generated;
    
	if (debug)
	  std::cerr << "generation: earley" << std::endl;
    
	utils::resource start;
    
	cicada::generate_earley(hypergraph, generated, depth, width, context);
    
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
      
      
      int depth;
      int width;
      bool context;
      
      int debug;
    };

  };
};


#endif
