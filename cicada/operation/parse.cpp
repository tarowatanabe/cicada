//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/parse.hpp>
#include <cicada/grammar_simple.hpp>

#include <cicada/operation/functional.hpp>
#include <cicada/operation/parse.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {

    ParseCKY::ParseCKY(const std::string& parameter,
		       const grammar_type& __grammar,
		       const std::string& __goal,
		       const std::string& __non_terminal,
		       const bool __insertion,
		       const bool __deletion,
		       const int __debug)
      : grammar(__grammar),
	goal(__goal), non_terminal(__non_terminal), 
	insertion(__insertion), deletion(__deletion),
	weights(0),
	size(200),
	weights_one(false),
	yield_source(false),
	treebank(false),
	debug(__debug)
    { 
      typedef cicada::Parameter param_type;
      
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "parse-cky" && utils::ipiece(param.name()) != "parse-cyk")
	throw std::runtime_error("this is not a CKY(CYK) parser");
      
      bool source = false;
      bool target = false;
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "size")
	  size = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "weights")
	  weights = &base_type::weights(piter->second);
	else if (utils::ipiece(piter->first) == "weights-one")
	  weights_one = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "yield") {
	  if (utils::ipiece(piter->second) == "source")
	    source = true;
	  else if (utils::ipiece(piter->second) == "target")
	    target = true;
	  else
	    throw std::runtime_error("unknown yield: " + piter->second);
	} else if (utils::ipiece(piter->first) == "treebank")
	  treebank = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for CKY parser: " << piter->first << "=" << piter->second << std::endl;
      }
	
      if (source && target)
	throw std::runtime_error("CKY parser can work either source or target yield");
	
      yield_source = source;

      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
    }

    void ParseCKY::assign(const weight_set_type& __weights)
    {
      weights = &__weights;
    }

    void ParseCKY::operator()(data_type& data) const
    {
      typedef cicada::semiring::Logprob<double> weight_type;

      const lattice_type& lattice = data.lattice;
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type parsed;
      
      hypergraph.clear();
      if (lattice.empty()) return;
    
      if (debug)
	std::cerr << "parse cky: " << data.id << std::endl;

      weight_set_type weights_zero;
      const weight_set_type* weights_parse = (weights ? weights : &weights_zero);

      utils::resource start;

      grammar_type grammar_parse(grammar);
    
      if (insertion)
	grammar_parse.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(lattice, non_terminal)));
      if (deletion)
	grammar_parse.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(lattice, non_terminal)));
	
      if (weights_one)
	cicada::parse_cky(goal, grammar_parse, weight_function_one<weight_type>(), lattice, parsed, size, yield_source, treebank);
      else
	cicada::parse_cky(goal, grammar_parse, weight_function<weight_type>(*weights_parse), lattice, parsed, size, yield_source, treebank);
      
      utils::resource end;
    
      if (debug)
	std::cerr << "parse cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << "parse: " << data.id
		  << " # of nodes: " << parsed.nodes.size()
		  << " # of edges: " << parsed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(parsed.is_valid())
		  << std::endl;
    
      hypergraph.swap(parsed);
    }
  };
};
