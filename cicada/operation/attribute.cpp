//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/push_head.hpp>

#include <cicada/operation/attribute.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    Attribute::Attribute(const std::string& parameter, const int __debug)
      :  base_type("attribute"), head_node(false), attr_head_node("head-node"), debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "attribute")
	throw std::runtime_error("this is not a head pusher");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	std::cerr << "WARNING: unsupported parameter for attribute: " << piter->first << "=" << piter->second << std::endl;
    }

    void Attribute::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type annotated;
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      utils::resource start;
      
      annotated = hypergraph;
      
      hypergraph_type::edge_set_type::iterator eiter_end = annotated.edges.end();
      for (hypergraph_type::edge_set_type::iterator eiter = annotated.edges.begin(); eiter != eiter_end; ++ eiter) {
	hypergraph_type::edge_type& edge = *eiter;
	
	edge.attributes[attr_head_node] = attribute_set_type::int_type(edge.head);
      }
      
      utils::resource end;
	
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << " thread time: " << (end.thread_time() - start.thread_time())
		  << std::endl;
	
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << annotated.nodes.size()
		  << " # of edges: " << annotated.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(annotated.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += annotated.nodes.size();
      stat.edge += annotated.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
      
      hypergraph.swap(annotated);
    }
  };
};
