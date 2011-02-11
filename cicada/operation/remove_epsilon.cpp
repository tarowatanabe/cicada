//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/parameter.hpp>
#include <cicada/remove_epsilon.hpp>

#include <cicada/operation/remove_epsilon.hpp>

#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    RemoveEpsilon::RemoveEpsilon(const std::string& parameter, const int __debug)
      :  debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "remove-epsilon")
	throw std::runtime_error("this is not an epsilon remover");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	std::cerr << "WARNING: unsupported parameter for remove-epsilon: " << piter->first << "=" << piter->second << std::endl;
    }

    void RemoveEpsilon::operator()(data_type& data) const
    {
      if (data.lattice.empty()) return;
	
      lattice_type removed;
	
      if (debug)
	std::cerr << "remove epsilon" << std::endl;
	
      utils::resource start;
	
      cicada::remove_epsilon(data.lattice, removed);
	
      utils::resource end;
	
      if (debug)
	std::cerr << "remove epsilon cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
	
      if (debug)
	std::cerr << "# of nodes: " << removed.size()
		  << " shortest distance: " << removed.shortest_distance()
		  << " longest distance: " << removed.longest_distance()
		  << std::endl;
	
      data.lattice.swap(removed);
    }
    
  };
};
