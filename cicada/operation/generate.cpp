//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/generate.hpp>

#include <cicada/operation/generate.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    GenerateEarley::GenerateEarley(const std::string& parameter,
				   const grammar_type& __grammar,
				   const std::string& __goal,
				   const std::string& __non_terminal,
				   const bool __insertion,
				   const bool __deletion,
				   const int __debug)
      : depth(0), width(0), debug(__debug)
    { 
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "generate-earley")
	throw std::runtime_error("this is not a Earley generator");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "depth")
	  depth = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "width")
	  width = utils::lexical_cast<int>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for generator: " << piter->first << "=" << piter->second << std::endl;
      }
    }

    void GenerateEarley::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;

      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type generated;
    
      if (debug)
	std::cerr << "generation: earley" << std::endl;
    
      utils::resource start;
    
      cicada::generate_earley(hypergraph, generated, depth, width);
    
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

  };
};
