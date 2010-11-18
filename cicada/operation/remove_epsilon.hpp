// -*- mode: c++ -*-

#ifndef __CICADA__OPERATION__REMOVE_EPSILON__HPP__
#define __CICADA__OPERATION__REMOVE_EPSILON__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/remove_epsilon.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

namespace cicada
{
  namespace operation
  {

    class RemoveEpsilon : public cicada::Operation
    {
    public:
      RemoveEpsilon(const std::string& parameter, const int __debug)
	:  debug(__debug)
      {
	typedef cicada::Parameter param_type;
	
	param_type param(parameter);
	if (param.name() != "remove-epsilon")
	  throw std::runtime_error("this is not an epsilon remover");
	
	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	  std::cerr << "WARNING: unsupported parameter for remove-epsilon: " << piter->first << "=" << piter->second << std::endl;
      }
  
  
      void operator()(data_type& data) const
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
  
      int debug;
    };

  };
};

#endif
