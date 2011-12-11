//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/remove_epsilon.hpp>

#include <cicada/operation/remove_epsilon.hpp>

#include <utils/resource.hpp>
#include <utils/piece.hpp>
#include <utils/lexical_cast.hpp>

namespace cicada
{
  namespace operation
  {
    RemoveEpsilon::RemoveEpsilon(const std::string& parameter, const int __debug)
      :  lattice_mode(false), forest_mode(false), debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "remove-epsilon")
	throw std::runtime_error("this is not an epsilon remover");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "lattice")
	  lattice_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "forest" || utils::ipiece(piter->first) == "hypergraph")
	  forest_mode = utils::lexical_cast<bool>(piter->second);
	else
	std::cerr << "WARNING: unsupported parameter for remove-epsilon: " << piter->first << "=" << piter->second << std::endl;
      }

      if (int(lattice_mode) + forest_mode > 1)
	throw std::runtime_error("specify either lattice or forest");
      if (int(lattice_mode) + forest_mode == 0)
	throw std::runtime_error("specify either lattice or forest");

      name = std::string("remove-epsilon-") + (lattice_mode ? "lattice" : "forest");
    }

    void RemoveEpsilon::operator()(data_type& data) const
    {
      if (lattice_mode) {
	if (data.lattice.empty()) return;
	
	lattice_type removed;
	
	if (debug)
	  std::cerr << name << ": " << data.id << std::endl;
	
	utils::resource start;
	
	cicada::remove_epsilon(data.lattice, removed);
	
	utils::resource end;
	
	if (debug)
	  std::cerr << name << ": " << data.id
		    << " cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << " thread time: " << (end.thread_time() - start.thread_time())
		    << std::endl;
	
	if (debug)
	  std::cerr << name << ": " << data.id
		    << " # of nodes: " << removed.size()
		    << " shortest distance: " << removed.shortest_distance()
		    << " longest distance: " << removed.longest_distance()
		    << std::endl;

	statistics_type::statistic_type& stat = data.statistics[name];
	
	++ stat.count;
	stat.node += removed.node_size();
	stat.edge += removed.edge_size();
	stat.user_time += (end.user_time() - start.user_time());
	stat.cpu_time  += (end.cpu_time() - start.cpu_time());
	stat.thread_time  += (end.thread_time() - start.thread_time());
	
	data.lattice.swap(removed);
      } else if (forest_mode) {
	if (! data.hypergraph.is_valid()) return;
	
	hypergraph_type removed;
	
	if (debug)
	  std::cerr << name << ": " << data.id << std::endl;
	
	utils::resource start;
	
	cicada::remove_epsilon(data.hypergraph, removed);
	
	utils::resource end;
	
	if (debug)
	  std::cerr << name << ": " << data.id
		    << " cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << " thread time: " << (end.thread_time() - start.thread_time())
		    << std::endl;
	
	if (debug)
	  std::cerr << name << ": " << data.id
		    << " # of nodes: " << removed.nodes.size()
		    << " # of edges: " << removed.edges.size()
		    << " valid? " << utils::lexical_cast<std::string>(removed.is_valid())
		    << std::endl;

	statistics_type::statistic_type& stat = data.statistics[name];
	
	++ stat.count;
	stat.node += removed.nodes.size();
	stat.edge += removed.edges.size();
	stat.user_time += (end.user_time() - start.user_time());
	stat.cpu_time  += (end.cpu_time() - start.cpu_time());
	stat.thread_time  += (end.thread_time() - start.thread_time());
	
	data.hypergraph.swap(removed);
	
      } else
	throw std::runtime_error("which one? lattice or forest");
    }
    
  };
};
