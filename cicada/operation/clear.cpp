//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>

#include <cicada/operation/clear.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    Clear::Clear(const std::string& parameter,
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
      if (utils::ipiece(param.name()) != "clear")
	throw std::runtime_error("this is not data clearer");
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "hypergraph")
	  clear_hypergraph = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "lattice")
	  clear_lattice = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "spans" || utils::ipiece(piter->first) == "span")
	  clear_spans = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "targets" || utils::ipiece(piter->first) == "bitext")
	  clear_targets = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "counts" || utils::ipiece(piter->first) == "ngram-counts")
	  clear_counts = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for clear: " << piter->first << "=" << piter->second << std::endl;
      }
    }
    
    void Clear::operator()(data_type& data) const
    {
      if (debug)
	std::cerr << "clear: " << data.id << std::endl;

      if (clear_hypergraph)
	data.hypergraph.clear();
      if (clear_lattice)
	data.lattice.clear();
      if (clear_spans)
	data.spans.clear();
      if (clear_targets)
	data.targets.clear();
      if (clear_counts)
	data.ngram_counts.clear();
    }
  };
};
