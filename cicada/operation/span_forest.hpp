// -*- mode: c++ -*-

#ifndef __CICADA__OPERATION__SPAN_FOREST__HPP__
#define __CICADA__OPERATION__SPAN_FOREST__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/span_forest.hpp>

#include <cicada/operation/functional.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

namespace cicada
{
  namespace operation
  {
    class SpanForest : public Operation
    {
    public:
      SpanForest(const std::string& parameter, const int __debug)
	: debug(__debug)
      {
	typedef cicada::Parameter param_type;
    
	param_type param(parameter);
	if (param.name() != "span-forest")
	  throw std::runtime_error("this is not a span annotator");

	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	  std::cerr << "WARNING: unsupported parameter for span annotator: " << piter->first << "=" << piter->second << std::endl;
      }

      void operator()(data_type& data) const
      {
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type spanned;

	if (debug)
	  std::cerr << "span annotation:"
		    << " # of nodes: " << hypergraph.nodes.size()
		    << " # of edges: " << hypergraph.edges.size()
		    << " valid? " << utils::lexical_cast<std::string>(hypergraph.is_valid())
		    << std::endl;
    
	utils::resource start;

	cicada::span_forest(hypergraph, spanned);
	
	utils::resource end;
    
	if (debug)
	  std::cerr << "span annotation cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << std::endl;
    
	if (debug)
	  std::cerr << "spanned:"
		    << " # of nodes: " << spanned.nodes.size()
		    << " # of edges: " << spanned.edges.size()
		    << " valid? " << utils::lexical_cast<std::string>(spanned.is_valid())
		    << std::endl;
    
	hypergraph.swap(spanned);

      }

      int debug;
    };

  };
};


#endif
