//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/span_forest.hpp>

#include <cicada/operation/span_forest.hpp>
#include <cicada/operation/functional.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    SpanForest::SpanForest(const std::string& parameter, const int __debug)
      : base_type("span-forest"), debug(__debug)
    {
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "span-forest")
	throw std::runtime_error("this is not a span annotator");

      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	std::cerr << "WARNING: unsupported parameter for span annotator: " << piter->first << "=" << piter->second << std::endl;
    }
    
    void SpanForest::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type spanned;

      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
    
      utils::resource start;

      cicada::span_forest(hypergraph, spanned);
	
      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << spanned.nodes.size()
		  << " # of edges: " << spanned.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(spanned.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
	
	++ stat.count;
	stat.node += spanned.nodes.size();
	stat.edge += spanned.edges.size();
	stat.user_time += (end.user_time() - start.user_time());
	stat.cpu_time  += (end.cpu_time() - start.cpu_time());
	stat.thread_time  += (end.thread_time() - start.thread_time());
    
      hypergraph.swap(spanned);
    }
  };
};
