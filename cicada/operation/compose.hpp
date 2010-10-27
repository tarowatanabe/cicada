// -*- mode: c++ -*-

#ifndef __CICADA__OPERATION__COMPOSE__HPP__
#define __CICADA__OPERATION__COMPOSE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/compose.hpp>
#include <cicada/grammar_hiero.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

namespace cicada
{
  namespace operation
  {
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
	  debug(__debug)
      {
	typedef cicada::Parameter param_type;
	
	param_type param(parameter);
	if (param.name() != "compose-earley")
	  throw std::runtime_error("this is not a Earley composer");
      }
  
      void operator()(data_type& data) const
      {
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type composed;
    
	if (debug)
	  std::cerr << "composition: earley" << std::endl;

	utils::resource start;

	grammar_type grammar_translation(grammar);
    
	if (insertion)
	  grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(hypergraph, non_terminal)));
	if (deletion)
	  grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(hypergraph, non_terminal)));

    
	cicada::compose_earley(grammar_translation, hypergraph, composed);
    
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
  
      const std::string goal;
      const std::string non_terminal;
  
      const bool insertion;
      const bool deletion;
  
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
	  debug(__debug)
      { 
	typedef cicada::Parameter param_type;
	
	param_type param(parameter);
	if (param.name() != "compose-cky")
	  throw std::runtime_error("this is not a CKY composer");
      }
  
      void operator()(data_type& data) const
      {
	const lattice_type& lattice = data.lattice;
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type composed;
    
	if (debug)
	  std::cerr << "composition: cky" << std::endl;

	utils::resource start;

	grammar_type grammar_translation(grammar);
    
	if (insertion)
	  grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(lattice, non_terminal)));
	if (deletion)
	  grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(lattice, non_terminal)));

    
	cicada::compose_cky(goal, grammar_translation, lattice, composed);
    
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
  
      const std::string goal;
      const std::string non_terminal;
  
      const bool insertion;
      const bool deletion;
  
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
	  debug(__debug)
      { 
	typedef cicada::Parameter param_type;
    
	param_type param(parameter);
	if (param.name() != "compose-phrase")
	  throw std::runtime_error("this is not a phrase composer");
	
	for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "distortion") == 0)
	    distortion = boost::lexical_cast<int>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for composer: " << piter->first << "=" << piter->second << std::endl;
	}
	
      }
  
      void operator()(data_type& data) const
      {
	const lattice_type& lattice = data.lattice;
	hypergraph_type& hypergraph = data.hypergraph;
	hypergraph_type composed;
    
	if (debug)
	  std::cerr << "composition: phrase" << std::endl;

	utils::resource start;

	grammar_type grammar_translation(grammar);
    
	if (insertion)
	  grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(lattice, non_terminal)));
	if (deletion)
	  grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(lattice, non_terminal)));

    
	cicada::compose_phrase(non_terminal, grammar_translation, lattice, distortion, composed);
    
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
  
      const std::string goal;
      const std::string non_terminal;
  
      const bool insertion;
      const bool deletion;

      int distortion;
  
      int debug;
    };
    
  };
};


#endif
