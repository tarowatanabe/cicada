//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/intersect.hpp>

#include <cicada/operation/intersect.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    Intersect::Intersect(const std::string& parameter, const int __debug)
      : lattice_mode(false), target_mode(false), debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "intersect")
	throw std::runtime_error("this is not a intersector");

      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "lattice")
	  lattice_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "target")
	  target_mode = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for intersect: " << piter->first << "=" << piter->second << std::endl;
      }

      if (lattice_mode && target_mode)
	throw std::runtime_error("either lattice or target");
	
      if (! lattice_mode && ! target_mode)
	target_mode = true;
      
      name = std::string("intersect-") + (lattice_mode ? "lattice" : "forest");

    }

    void Intersect::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;
      
      if (lattice_mode) {
	if (data.lattice.empty()) {
	  data.hypergraph.clear();
	  return;
	}
      } else {
	if (data.targets.empty() || data.targets.front().empty()) {
	  data.hypergraph.clear();
	  return;
	}
      }
      
      const sentence_set_type& targets = data.targets;
      hypergraph_type& hypergraph = data.hypergraph;
      
      hypergraph_type intersected;

      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      utils::resource start;
      
      if (lattice_mode)
	cicada::intersect(hypergraph, data.lattice, intersected);
      else {
	lattice_type target(targets.front());
	
	cicada::intersect(hypergraph, target, intersected);
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
		  << " # of nodes: " << intersected.nodes.size()
		  << " # of edges: " << intersected.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(intersected.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += intersected.nodes.size();
      stat.edge += intersected.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
      
      hypergraph.swap(intersected);
    }

  };
};
