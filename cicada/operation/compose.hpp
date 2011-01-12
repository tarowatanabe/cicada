// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__COMPOSE__HPP__
#define __CICADA__OPERATION__COMPOSE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/compose.hpp>
#include <cicada/grammar_simple.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

namespace cicada
{
  namespace operation
  {
    class ComposeTree : public Operation
    {
    public:
      ComposeTree(const std::string& parameter,
		  const tree_grammar_type& __tree_grammar,
		  const grammar_type& __grammar,
		  const std::string& __goal,
		  const std::string& __non_terminal,
		  const bool __insertion,
		  const bool __deletion,
		  const bool __fallback,
		  const int __debug)
	: tree_grammar(__tree_grammar), grammar(__grammar),
	  goal(__goal), non_terminal(__non_terminal), 
	  insertion(__insertion), deletion(__deletion), fallback(__fallback),
	  yield_source(false),
	  debug(__debug)
      {
	typedef cicada::Parameter param_type;
	
	param_type param(parameter);
	if (param.name() != "compose-tree")
	  throw std::runtime_error("this is not a Tree composer");
	
	bool source = false;
	bool target = false;
	
	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	    if (strcasecmp(piter->second.c_str(), "source") == 0)
	      source = true;
	    else if (strcasecmp(piter->second.c_str(), "target") == 0)
	      target = true;
	    else
	      throw std::runtime_error("unknown yield: " + piter->second);
	  } else
	    std::cerr << "WARNING: unsupported parameter for Tree composer: " << piter->first << "=" << piter->second << std::endl;
	}
	
	if (source && target)
	  throw std::runtime_error("Tree composer can work either source or target yield");
	
	yield_source = source;
      }
      
      void operator()(data_type& data) const
      {
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type composed;
    
	if (debug)
	  std::cerr << "composition: tree" << std::endl;
	
	utils::resource start;
	
	grammar_type grammar_compose(grammar);
    
	if (insertion)
	  grammar_compose.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(hypergraph, non_terminal)));
	if (deletion)
	  grammar_compose.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(hypergraph, non_terminal)));


	tree_grammar_type tree_grammar_compose(tree_grammar);
	if (fallback)
	  tree_grammar_compose.push_back(tree_grammar_type::transducer_ptr_type(new cicada::TreeGrammarFallback(hypergraph, non_terminal)));
	
	cicada::compose_tree(goal, tree_grammar_compose, grammar_compose, hypergraph, composed, yield_source);
	
	utils::resource end;
    
	if (debug)
	  std::cerr << "compose cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << std::endl;
    
	if (debug)
	  std::cerr << "# of nodes: " << composed.nodes.size()
		    << " # of edges: " << composed.edges.size()
		    << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		    << std::endl;
	
	hypergraph.swap(composed);
      }
      
      const tree_grammar_type& tree_grammar;
      const grammar_type&      grammar;
      
      std::string goal;
      std::string non_terminal;
      
      bool insertion;
      bool deletion;
      bool fallback;
      
      bool yield_source;
      
      int debug;
    };

    class ComposeEarley : public Operation
    {
    public:
      ComposeEarley(const std::string& parameter,
		    const grammar_type& __grammar,
		    const std::string& __goal,
		    const std::string& __non_terminal,
		    const bool __insertion,
		    const bool __deletion,
		    const int __debug)
	: grammar(__grammar),
	  goal(__goal), non_terminal(__non_terminal), 
	  insertion(__insertion), deletion(__deletion),
	  yield_source(false),
	  debug(__debug)
      {
	typedef cicada::Parameter param_type;
	
	param_type param(parameter);
	if (param.name() != "compose-earley")
	  throw std::runtime_error("this is not a Earley composer");
	
	bool source = false;
	bool target = false;
	
	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	    if (strcasecmp(piter->second.c_str(), "source") == 0)
	      source = true;
	    else if (strcasecmp(piter->second.c_str(), "target") == 0)
	      target = true;
	    else
	      throw std::runtime_error("unknown yield: " + piter->second);
	  } else
	    std::cerr << "WARNING: unsupported parameter for Earley composer: " << piter->first << "=" << piter->second << std::endl;
	}
	
	if (source && target)
	  throw std::runtime_error("Earley composer can work either source or target yield");
	
	yield_source = source;
      }
  
      void operator()(data_type& data) const
      {
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type composed;
    
	if (debug)
	  std::cerr << "composition: earley" << std::endl;

	utils::resource start;

	grammar_type grammar_compose(grammar);
    
	if (insertion)
	  grammar_compose.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(hypergraph, non_terminal)));
	if (deletion)
	  grammar_compose.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(hypergraph, non_terminal)));

    
	cicada::compose_earley(grammar_compose, hypergraph, composed, yield_source);
    
	utils::resource end;
    
	if (debug)
	  std::cerr << "compose cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << std::endl;
    
	if (debug)
	  std::cerr << "# of nodes: " << composed.nodes.size()
		    << " # of edges: " << composed.edges.size()
		    << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		    << std::endl;
    
	hypergraph.swap(composed);
      }
  
      const grammar_type& grammar;
  
      std::string goal;
      std::string non_terminal;
  
      bool insertion;
      bool deletion;

      bool yield_source;
  
      int debug;
    };

    class ComposeCKY : public Operation
    {
    public:
      ComposeCKY(const std::string& parameter,
		 const grammar_type& __grammar,
		 const std::string& __goal,
		 const std::string& __non_terminal,
		 const bool __insertion,
		 const bool __deletion,
		 const int __debug)
	: grammar(__grammar),
	  goal(__goal), non_terminal(__non_terminal), 
	  insertion(__insertion), deletion(__deletion),
	  yield_source(false),
	  debug(__debug)
      { 
	typedef cicada::Parameter param_type;
	
	param_type param(parameter);
	if (param.name() != "compose-cky")
	  throw std::runtime_error("this is not a CKY composer");

	bool source = false;
	bool target = false;
	
	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	    if (strcasecmp(piter->second.c_str(), "source") == 0)
	      source = true;
	    else if (strcasecmp(piter->second.c_str(), "target") == 0)
	      target = true;
	    else
	      throw std::runtime_error("unknown yield: " + piter->second);
	  } else
	    std::cerr << "WARNING: unsupported parameter for CKY composer: " << piter->first << "=" << piter->second << std::endl;
	}
	
	if (source && target)
	  throw std::runtime_error("CKY composer can work either source or target yield");
	
	yield_source = source;
      }
  
      void operator()(data_type& data) const
      {
	const lattice_type& lattice = data.lattice;
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type composed;
    
	if (debug)
	  std::cerr << "composition: cky" << std::endl;

	utils::resource start;

	grammar_type grammar_compose(grammar);
    
	if (insertion)
	  grammar_compose.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(lattice, non_terminal)));
	if (deletion)
	  grammar_compose.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(lattice, non_terminal)));
	
	cicada::compose_cky(goal, grammar_compose, lattice, composed, yield_source);
    
	utils::resource end;
    
	if (debug)
	  std::cerr << "compose cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << std::endl;
    
	if (debug)
	  std::cerr << "# of nodes: " << composed.nodes.size()
		    << " # of edges: " << composed.edges.size()
		    << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		    << std::endl;
    
	hypergraph.swap(composed);
      }
  
      const grammar_type& grammar;
  
      std::string goal;
      std::string non_terminal;
  
      bool insertion;
      bool deletion;

      bool yield_source;
  
      int debug;
    };


    class ComposePhrase : public Operation
    {
    public:
      ComposePhrase(const std::string& parameter,
		    const grammar_type& __grammar,
		    const std::string& __goal,
		    const std::string& __non_terminal,
		    const bool __insertion,
		    const bool __deletion,
		    const int __debug)
	: grammar(__grammar),
	  goal(__goal), non_terminal(__non_terminal), 
	  insertion(__insertion), deletion(__deletion),
	  distortion(0),
	  yield_source(false),
	  debug(__debug)
      { 
	typedef cicada::Parameter param_type;
    
	param_type param(parameter);
	if (param.name() != "compose-phrase")
	  throw std::runtime_error("this is not a phrase composer");

	bool source = false;
	bool target = false;
	
	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "distortion") == 0)
	    distortion = boost::lexical_cast<int>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	    if (strcasecmp(piter->second.c_str(), "source") == 0)
	      source = true;
	    else if (strcasecmp(piter->second.c_str(), "target") == 0)
	      target = true;
	    else
	      throw std::runtime_error("unknown yield: " + piter->second);
	  } else
	    std::cerr << "WARNING: unsupported parameter for composer: " << piter->first << "=" << piter->second << std::endl;
	}
	
	if (source && target)
	  throw std::runtime_error("Phrase composer can work either source or target yield");
	
	yield_source = source;
      }
  
      void operator()(data_type& data) const
      {
	const lattice_type& lattice = data.lattice;
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type composed;
    
	if (debug)
	  std::cerr << "composition: phrase" << std::endl;

	utils::resource start;

	grammar_type grammar_compose(grammar);
    
	if (insertion)
	  grammar_compose.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(lattice, non_terminal)));
	if (deletion)
	  grammar_compose.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(lattice, non_terminal)));

    
	cicada::compose_phrase(non_terminal, grammar_compose, lattice, distortion, composed);
    
	utils::resource end;
    
	if (debug)
	  std::cerr << "compose cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << std::endl;
    
	if (debug)
	  std::cerr << "# of nodes: " << composed.nodes.size()
		    << " # of edges: " << composed.edges.size()
		    << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		    << std::endl;
    
	hypergraph.swap(composed);
      }
  
      const grammar_type& grammar;
  
      std::string goal;
      std::string non_terminal;
  
      bool insertion;
      bool deletion;

      int distortion;
      
      bool yield_source;
  
      int debug;
    };


    class ComposeAlignment : public Operation
    {
    public:
      ComposeAlignment(const std::string& parameter,
		       const grammar_type& __grammar,
		       const std::string& __goal,
		       const std::string& __non_terminal,
		       const int __debug)
	: grammar(__grammar),
	  goal(__goal), non_terminal(__non_terminal), 
	  lattice_mode(false),
	  forest_mode(false),
	  debug(__debug)
      { 
	typedef cicada::Parameter param_type;
    
	param_type param(parameter);
	if (param.name() != "compose-alignment")
	  throw std::runtime_error("this is not a alignment composer");

	bool source = false;
	bool target = false;
	
	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "lattice") == 0)
	    lattice_mode = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "forest") == 0) {
	    forest_mode = utils::lexical_cast<bool>(piter->second);
	  } else
	    std::cerr << "WARNING: unsupported parameter for composer: " << piter->first << "=" << piter->second << std::endl;
	}
	
	if (lattice_mode && forest_mode)
	  throw std::runtime_error("either lattice or forest");

	if (! lattice_mode && ! forest_mode)
	  lattice_mode = true;
      }
  
      void operator()(data_type& data) const
      {
	const lattice_type& lattice = data.lattice;
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type composed;
	
	lattice_type target;
	if (! data.targets.empty())
	  target = lattice_type(data.targets.front());
	
	if (debug)
	  std::cerr << "composition: alignment" << std::endl;
	
	utils::resource start;
	
	grammar_type grammar_alignment(grammar);
	if (lattice_mode)
	  grammar_alignment.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarPair(lattice, target, non_terminal)));
	else
	  grammar_alignment.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarPair(hypergraph, target, non_terminal)));
	
	if (lattice_mode)
	  cicada::compose_alignment(non_terminal, grammar_alignment, lattice, target, composed);
	else
	  cicada::compose_alignment(non_terminal, grammar_alignment, hypergraph, target, composed);
    
	utils::resource end;
    
	if (debug)
	  std::cerr << "compose cpu time: " << (end.cpu_time() - start.cpu_time())
		    << " user time: " << (end.user_time() - start.user_time())
		    << std::endl;
    
	if (debug)
	  std::cerr << "# of nodes: " << composed.nodes.size()
		    << " # of edges: " << composed.edges.size()
		    << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		    << std::endl;
    
	hypergraph.swap(composed);
      }
  
      const grammar_type& grammar;
  
      std::string goal;
      std::string non_terminal;

      bool lattice_mode;
      bool forest_mode;
  
      int debug;
    };
    
  };
};


#endif
