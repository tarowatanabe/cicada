// -*- mode: c++ -*-

#ifndef __CICADA__OPERATION__CLEAR__HPP__
#define __CICADA__OPERATION__CLEAR__HPP__ 1

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
    class Clear : public Operation
    {
    public:
      Clear(const std::string& parameter,
	    const int __debug)
	: clear_hypergraph(false),
	  clear_lattice(false),
	  clear_spans(false),
	  clear_targets(false),
	  clear_counts(false),
	  debug(__debug)
      { 
	typedef cicada::Parameter param_type;
    
	param_type param(parameter);
	if (param.name() != "generate-earley")
	  throw std::runtime_error("this is not a Earley generator");
	
	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "hypergraph") == 0)
	    clear_hypergraph = utils::lexical_cast<bool>(piter->second);
	  if (strcasecmp(piter->first.c_str(), "lattice") == 0)
	    clear_lattice = utils::lexical_cast<bool>(piter->second);
	  if (strcasecmp(piter->first.c_str(), "spans") == 0 || strcasecmp(piter->first.c_str(), "span") == 0)
	    clear_spans = utils::lexical_cast<bool>(piter->second);
	  if (strcasecmp(piter->first.c_str(), "targets") == 0 || strcasecmp(piter->first.c_str(), "bitext") == 0)
	    clear_targets = utils::lexical_cast<bool>(piter->second);
	  if (strcasecmp(piter->first.c_str(), "counts") == 0 || strcasecmp(piter->first.c_str(), "ngram-counts") == 0)
	    clear_counts = utils::lexical_cast<bool>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for generator: " << piter->first << "=" << piter->second << std::endl;
	}
      }
  
      void operator()(data_type& data) const
      {
	if (debug)
	  std::cerr << "clear" << std::endl;
	
      }
      
      bool clear_hypergraph;
      bool clear_lattice;
      bool clear_spans;
      bool clear_targets;
      bool clear_counts;
      
      int debug;
    };

  };
};


#endif
