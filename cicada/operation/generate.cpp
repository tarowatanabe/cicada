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
				   const int __debug)
      : base_type("generate-earley"),
	depth(0), width(0), debug(__debug)
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
	std::cerr << name << ": " << data.id << std::endl;
      
      utils::resource start;
      
      cicada::generate_earley(hypergraph, generated, depth, width);
      
      utils::resource end;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << generated.nodes.size()
		  << " # of edges: " << generated.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(generated.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += generated.nodes.size();
      stat.edge += generated.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
    
      hypergraph.swap(generated);
    }

  };
};
