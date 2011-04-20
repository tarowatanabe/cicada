//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/parse.hpp>
#include <cicada/grammar_unknown.hpp>

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
		       const int __debug)
      : base_type("parse-cky"),
	grammar(__grammar),
	goal(__goal), 
	weights(0),
	weights_assigned(0),
	size(200),
	weights_one(false),
	weights_fixed(false),
	yield_source(false),
	treebank(false),
	pos_mode(false),
	unique_goal(false),
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
	else if (utils::ipiece(piter->first) == "goal")
	  goal = piter->second;
	else if (utils::ipiece(piter->first) == "grammar")
	  grammar_local.push_back(piter->second);
	else if (utils::ipiece(piter->first) == "yield") {
	  if (utils::ipiece(piter->second) == "source")
	    source = true;
	  else if (utils::ipiece(piter->second) == "target")
	    target = true;
	  else
	    throw std::runtime_error("unknown yield: " + piter->second);
	} else if (utils::ipiece(piter->first) == "treebank")
	  treebank = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "pos")
	  pos_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "unique" || utils::ipiece(piter->first) == "unique-goal")
	  unique_goal = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for CKY parser: " << piter->first << "=" << piter->second << std::endl;
      }
	
      if (source && target)
	throw std::runtime_error("CKY parser can work either source or target yield");
	
      yield_source = source;

      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
      
      if (weights || weights_one)
	weights_fixed = true;
      
      if (! weights)
	weights = &base_type::weights();
    }

    void ParseCKY::assign(const weight_set_type& __weights)
    {
      if (! weights_fixed)
	weights_assigned = &__weights;
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
	std::cerr << name << ": " << data.id << std::endl;

      const weight_set_type* weights_parse = (weights_assigned ? weights_assigned : &(weights->weights));
      
      const grammar_type& grammar_parse = (grammar_local.empty() ? grammar : grammar_local);

      utils::resource start;
      
      grammar_parse.assign(lattice);
	
      if (weights_one)
	cicada::parse_cky(goal, grammar_parse, weight_function_one<weight_type>(), lattice, parsed, size, yield_source, treebank, pos_mode, unique_goal);
      else
	cicada::parse_cky(goal, grammar_parse, weight_function<weight_type>(*weights_parse), lattice, parsed, size, yield_source, treebank, pos_mode, unique_goal);
      
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

    ParseAgenda::ParseAgenda(const std::string& parameter,
		       const grammar_type& __grammar,
		       const std::string& __goal,
		       const int __debug)
      : base_type("parse-agenda"),
	grammar(__grammar),
	goal(__goal),
	weights(0),
	weights_assigned(0),
	size(200),
	weights_one(false),
	weights_fixed(false),
	yield_source(false),
	treebank(false),
	pos_mode(false),
	debug(__debug)
    { 
      typedef cicada::Parameter param_type;
      
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "parse-agenda")
	throw std::runtime_error("this is not an agenda parser");
      
      bool source = false;
      bool target = false;
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "size")
	  size = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "weights")
	  weights = &base_type::weights(piter->second);
	else if (utils::ipiece(piter->first) == "weights-one")
	  weights_one = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "goal")
	  goal = piter->second;
	else if (utils::ipiece(piter->first) == "grammar")
	  grammar_local.push_back(piter->second);
	else if (utils::ipiece(piter->first) == "yield") {
	  if (utils::ipiece(piter->second) == "source")
	    source = true;
	  else if (utils::ipiece(piter->second) == "target")
	    target = true;
	  else
	    throw std::runtime_error("unknown yield: " + piter->second);
	} else if (utils::ipiece(piter->first) == "treebank")
	  treebank = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "pos")
	  pos_mode = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for agenda parser: " << piter->first << "=" << piter->second << std::endl;
      }
	
      if (source && target)
	throw std::runtime_error("agenda parser can work either source or target yield");
      
      yield_source = source;

      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
      
      if (weights || weights_one)
	weights_fixed = true;

      if (! weights)
	weights = &base_type::weights();
    }

    void ParseAgenda::assign(const weight_set_type& __weights)
    {
      if (! weights_fixed)
	weights_assigned = &__weights;
    }

    void ParseAgenda::operator()(data_type& data) const
    {
      typedef cicada::semiring::Logprob<double> weight_type;
      
      const lattice_type& lattice = data.lattice;
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type parsed;
      
      hypergraph.clear();
      if (lattice.empty()) return;
    
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      const weight_set_type* weights_parse = (weights_assigned ? weights_assigned : &(weights->weights));

      const grammar_type& grammar_parse = (grammar_local.empty() ? grammar : grammar_local);

      utils::resource start;

      grammar_parse.assign(lattice);
      
      if (weights_one)
	cicada::parse_agenda(goal, grammar_parse, weight_function_one<weight_type>(), lattice, parsed, size, yield_source, treebank, pos_mode);
      else
	cicada::parse_agenda(goal, grammar_parse, weight_function<weight_type>(*weights_parse), lattice, parsed, size, yield_source, treebank, pos_mode);
      
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
    
    
    ParseCoarse::ParseCoarse(const std::string& parameter,
			     const grammar_type& __grammar,
			     const std::string& __goal,
			     const int __debug)
      : base_type("parse-coarse"),
	grammars(), thresholds(), grammar(__grammar),
	goal(__goal),
	weights(0),
	weights_assigned(0),
	weights_one(false),
	weights_fixed(false),
	yield_source(false),
	treebank(false),
	pos_mode(false),
	debug(__debug)
    { 
      typedef cicada::Parameter param_type;
      
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "parse-coarse" && utils::ipiece(param.name()) != "parse-coarse")
	throw std::runtime_error("this is not an coarse parser");
      
      bool source = false;
      bool target = false;
      
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "weights")
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
	else if (utils::ipiece(piter->first) == "pos")
	  pos_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "goal")
	  goal = piter->second;
	else if (utils::ipiece(piter->first) == "grammar")
	  grammar_local.push_back(piter->second);
	else {
	  namespace qi = boost::spirit::qi;
	  
	  std::string::const_iterator iter = piter->first.begin();
	  std::string::const_iterator iter_end = piter->first.end();
	  
	  int id = -1;
	  if (qi::parse(iter, iter_end, "coarse" >> qi::int_, id) && iter == iter_end) {
	    if (id >= static_cast<int>(grammars.size()))
	      grammars.resize(id + 1);
	    
	    grammars[id].push_back(piter->second);
	  } else if (qi::parse(iter, iter_end, "threshold" >> qi::int_, id) && iter == iter_end) {
	    if (id >= static_cast<int>(thresholds.size()))
	      thresholds.resize(id + 1);
	    thresholds[id] = utils::lexical_cast<double>(piter->second);
	  } else
	    std::cerr << "WARNING: unsupported parameter for coarse parser: " << piter->first << "=" << piter->second << std::endl;
	}
      }
	
      if (source && target)
	throw std::runtime_error("coarse parser can work either source or target yield");
      
      yield_source = source;

      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
      
      if (weights || weights_one)
	weights_fixed = true;

      if (! weights)
	weights = &base_type::weights();
      
      if (grammars.size() != thresholds.size())
	throw std::runtime_error("# of coarse grammars and # of thresholds do not match");
      if (grammars.empty())
	throw std::runtime_error("no coarse grammar(s)?");
      
      // use either local or globally assigned grammar
      if (! grammar_local.empty())
	grammars.push_back(grammar_local);
      else
	grammars.push_back(grammar);
      
      // assign unknown grammar from the fine-grammar
      if (grammars.back().size() >= 2) {
	grammar_type::transducer_ptr_type unknown;
	
	grammar_type::iterator giter_end = grammars.back().end();
	for (grammar_type::iterator giter = grammars.back().begin(); giter != giter_end; ++ giter)
	  if (dynamic_cast<GrammarUnknown*>(&(*(*giter))))
	    unknown = *giter;
	
	if (unknown) {
	  grammar_set_type::iterator giter_end = grammars.end() - 1;
	  for (grammar_set_type::iterator giter = grammars.begin(); giter != giter_end; ++ giter)
	    if (giter->size() == 1)
	      giter->push_back(unknown);
	}
      }
    }

    void ParseCoarse::assign(const weight_set_type& __weights)
    {
      if (! weights_fixed)
	weights_assigned = &__weights;
    }

    void ParseCoarse::operator()(data_type& data) const
    {
      typedef cicada::semiring::Logprob<double> weight_type;
      
      const lattice_type& lattice = data.lattice;
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type parsed;
      
      hypergraph.clear();
      if (lattice.empty()) return;
    
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      const weight_set_type* weights_parse = (weights_assigned ? weights_assigned : &(weights->weights));

      utils::resource start;

      grammar_set_type::const_iterator giter_end = grammars.end();
      for (grammar_set_type::const_iterator giter = grammars.begin(); giter != giter_end; ++ giter)
	giter->assign(lattice);
      
      if (weights_one)
	cicada::parse_coarse(goal,
			     grammars.begin(), grammars.end(),
			     thresholds.begin(), thresholds.end(),
			     weight_function_one<weight_type>(), lattice, parsed, yield_source, treebank, pos_mode);
      else
	cicada::parse_coarse(goal,
			     grammars.begin(), grammars.end(),
			     thresholds.begin(), thresholds.end(),
			     weight_function<weight_type>(*weights_parse), lattice, parsed, yield_source, treebank, pos_mode);
      
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
