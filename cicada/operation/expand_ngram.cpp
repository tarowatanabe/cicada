//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/expand_ngram.hpp>

#include <cicada/operation/expand_ngram.hpp>
#include <cicada/operation/functional.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    ExpandNGram::ExpandNGram(const std::string& parameter, const int __debug)
      : base_type("expand-ngram"),
	order(0),
	debug(__debug)
    {
      typedef cicada::Parameter param_type;
      
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "expand-ngram")
	throw std::runtime_error("this is not an expand-ngram computer"); 
    
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (order <= 0)
	throw std::runtime_error("order must be positive");
    }

    void ExpandNGram::operator()(data_type& data) const
    {
      hypergraph_type& hypergraph = data.hypergraph;
      
      if (! hypergraph.is_valid()) return;
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      hypergraph_type expanded;
          
      utils::resource start;
      
      cicada::expand_ngram(hypergraph, expanded, order);

      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << " thread time: " << (end.thread_time() - start.thread_time())
		  << std::endl;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << expanded.nodes.size()
		  << " # of edges: " << expanded.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(expanded.is_valid())
		  << std::endl;
      
      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += expanded.nodes.size();
      stat.edge += expanded.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
      
      hypergraph.swap(expanded);
    }

  };
};
