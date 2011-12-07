//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/compose.hpp>

#include <cicada/operation/compose.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  namespace operation
  {
    ComposeTreeCKY::ComposeTreeCKY(const std::string& parameter,
				   const tree_grammar_type& __tree_grammar,
				   const grammar_type& __grammar,
				   const std::string& __goal,
				   const int __debug)
      : base_type("compose-tree-cky"),
	tree_grammar(__tree_grammar), grammar(__grammar),
	goal(__goal),
	yield_source(false),
	unique_goal(false), 
	debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "compose-tree-cky"
	  && utils::ipiece(param.name()) != "compose-tree-cyk"
	  && utils::ipiece(param.name()) != "compose-cky-tree"
	  && utils::ipiece(param.name()) != "compose-cyk-tree")
	throw std::runtime_error("this is not a Tree composer");
	
      bool source = false;
      bool target = false;
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "yield") {
	  if (utils::ipiece(piter->second) == "source")
	    source = true;
	  else if (utils::ipiece(piter->second) == "target")
	    target = true;
	  else
	    throw std::runtime_error("unknown yield: " + piter->second);
	} else if (utils::ipiece(piter->first) == "goal")
	  goal = piter->second;
	else if (utils::ipiece(piter->first) == "grammar")
	  grammar_local.push_back(piter->second);
	else if (utils::ipiece(piter->first) == "tree-grammar")
	  tree_grammar_local.push_back(piter->second);
	else if (utils::ipiece(piter->first) == "unique" || utils::ipiece(piter->first) == "unique-goal")
	  unique_goal = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for Tree composer: " << piter->first << "=" << piter->second << std::endl;
      }
	
      if (source && target)
	throw std::runtime_error("Tree composer can work either source or target yield");
	
      yield_source = source;
    }

    void ComposeTreeCKY::operator()(data_type& data) const
    {
      if (data.lattice.empty()) return;

      const lattice_type& lattice = data.lattice;
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type composed;

      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      const grammar_type& grammar_compose = (grammar_local.empty() ? grammar : grammar_local);
      const tree_grammar_type& tree_grammar_compose = (tree_grammar_local.empty() ? tree_grammar : tree_grammar_local);
	
      utils::resource start;

      grammar_compose.assign(lattice);
      tree_grammar_compose.assign(lattice);
      
      cicada::compose_tree_cky(goal, tree_grammar_compose, grammar_compose, lattice, composed, yield_source, unique_goal);
	
      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << composed.nodes.size()
		  << " # of edges: " << composed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += composed.nodes.size();
      stat.edge += composed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
      
      hypergraph.swap(composed);
    }

    ComposeTree::ComposeTree(const std::string& parameter,
			     const tree_grammar_type& __tree_grammar,
			     const grammar_type& __grammar,
			     const std::string& __goal,
			     const int __debug)
      : base_type("compose-tree"),
	tree_grammar(__tree_grammar), grammar(__grammar),
	goal(__goal),
	yield_source(false),
	debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "compose-tree")
	throw std::runtime_error("this is not a Tree composer");
	
      bool source = false;
      bool target = false;
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "yield") {
	  if (utils::ipiece(piter->second) == "source")
	    source = true;
	  else if (utils::ipiece(piter->second) == "target")
	    target = true;
	  else
	    throw std::runtime_error("unknown yield: " + piter->second);
	} else if (utils::ipiece(piter->first) == "goal")
	  goal = piter->second;
	else if (utils::ipiece(piter->first) == "grammar")
	  grammar_local.push_back(piter->second);
	else if (utils::ipiece(piter->first) == "tree-grammar")
	  tree_grammar_local.push_back(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for Tree composer: " << piter->first << "=" << piter->second << std::endl;
      }
	
      if (source && target)
	throw std::runtime_error("Tree composer can work either source or target yield");
	
      yield_source = source;
    }

    void ComposeTree::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;

      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type composed;

      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      const grammar_type& grammar_compose = (grammar_local.empty() ? grammar : grammar_local);
      const tree_grammar_type& tree_grammar_compose = (tree_grammar_local.empty() ? tree_grammar : tree_grammar_local);
	
      utils::resource start;

      grammar_compose.assign(hypergraph);
      tree_grammar_compose.assign(hypergraph);
	
      cicada::compose_tree(goal, tree_grammar_compose, grammar_compose, hypergraph, composed, yield_source);
	
      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << composed.nodes.size()
		  << " # of edges: " << composed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += composed.nodes.size();
      stat.edge += composed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
	
      hypergraph.swap(composed);
    }
    
    ComposeEarley::ComposeEarley(const std::string& parameter,
				 const grammar_type& __grammar,
				 const std::string& __goal,
				 const int __debug)
      : base_type("compose-earley"),
	grammar(__grammar),
	goal(__goal),
	yield_source(false),
	debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "compose-earley")
	throw std::runtime_error("this is not a Earley composer");
	
      bool source = false;
      bool target = false;
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "yield") {
	  if (utils::ipiece(piter->second) == "source")
	    source = true;
	  else if (utils::ipiece(piter->second) == "target")
	    target = true;
	  else
	    throw std::runtime_error("unknown yield: " + piter->second);
	} else if (utils::ipiece(piter->first) == "grammar")
	  grammar_local.push_back(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for Earley composer: " << piter->first << "=" << piter->second << std::endl;
      }
	
      if (source && target)
	throw std::runtime_error("Earley composer can work either source or target yield");
	
      yield_source = source;
    }
    
    void ComposeEarley::operator()(data_type& data) const
    {
      if (! data.hypergraph.is_valid()) return;

      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type composed;
    
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;

      const grammar_type& grammar_compose = (grammar_local.empty() ? grammar : grammar_local);

      utils::resource start;

      grammar_compose.assign(hypergraph);
    
      cicada::compose_earley(grammar_compose, hypergraph, composed, yield_source);
    
      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << composed.nodes.size()
		  << " # of edges: " << composed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += composed.nodes.size();
      stat.edge += composed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
    
      hypergraph.swap(composed);
    }

    ComposeCKY::ComposeCKY(const std::string& parameter,
			   const grammar_type& __grammar,
			   const std::string& __goal,
			   const int __debug)
      : base_type("compose-cky"),
	grammar(__grammar),
	goal(__goal),
	yield_source(false),
	treebank(false),
	pos_mode(false),
	unique_goal(false),
	debug(__debug)
    { 
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "compose-cky" && utils::ipiece(param.name()) != "compose-cyk")
	throw std::runtime_error("this is not a CKY(CYK) composer");

      bool source = false;
      bool target = false;
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "yield") {
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
	else if (utils::ipiece(piter->first) == "goal")
	  goal = piter->second;
	else if (utils::ipiece(piter->first) == "grammar")
	  grammar_local.push_back(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for CKY composer: " << piter->first << "=" << piter->second << std::endl;
      }
	
      if (source && target)
	throw std::runtime_error("CKY composer can work either source or target yield");
	
      yield_source = source;
    }

    void ComposeCKY::operator()(data_type& data) const
    {
      const lattice_type& lattice = data.lattice;
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type composed;
      
      hypergraph.clear();
      if (lattice.empty()) return;
    
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;

      const grammar_type& grammar_compose = (grammar_local.empty() ? grammar : grammar_local);

      utils::resource start;

      grammar_compose.assign(lattice);
	
      cicada::compose_cky(goal, grammar_compose, lattice, composed, yield_source, treebank, pos_mode, unique_goal);
    
      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << composed.nodes.size()
		  << " # of edges: " << composed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += composed.nodes.size();
      stat.edge += composed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
    
      hypergraph.swap(composed);
    }
    
    ComposeGrammar::ComposeGrammar(const std::string& parameter,
				   const grammar_type& __grammar,
				   const std::string& __goal,
				   const int __debug)
      : base_type("compose-grammar"),
	grammar(__grammar),
	goal(__goal),
	yield_source(false),
	debug(__debug)
    {
     
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "compose-grammar")
	throw std::runtime_error("this is not a grammar matching composer");

      bool source = false;
      bool target = false;
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "yield") {
	  if (utils::ipiece(piter->second) == "source")
	    source = true;
	  else if (utils::ipiece(piter->second) == "target")
	    target = true;
	  else
	    throw std::runtime_error("unknown yield: " + piter->second);
	} else if (utils::ipiece(piter->first) == "grammar")
	  grammar_local.push_back(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for composer: " << piter->first << "=" << piter->second << std::endl;
      }
	
      if (source && target)
	throw std::runtime_error("Phrase composer can work either source or target yield");
	
      yield_source = source;
    }

    void ComposeGrammar::operator()(data_type& data) const
    {
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type composed;

      if (! hypergraph.is_valid()) return;
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;

      const grammar_type& grammar_compose = (grammar_local.empty() ? grammar : grammar_local);

      utils::resource start;

      grammar_compose.assign(hypergraph);
      
      cicada::compose_grammar(grammar_compose, hypergraph, composed, yield_source);
    
      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << composed.nodes.size()
		  << " # of edges: " << composed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += composed.nodes.size();
      stat.edge += composed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
    
      hypergraph.swap(composed);
    }

    
    ComposePhrase::ComposePhrase(const std::string& parameter,
				 const grammar_type& __grammar,
				 const std::string& __goal,
				 const int __debug)
      : base_type("compose-phrase"),
	grammar(__grammar),
	goal(__goal),
	distortion(0),
	yield_source(false),
	debug(__debug)
    { 
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "compose-phrase")
	throw std::runtime_error("this is not a phrase composer");

      bool source = false;
      bool target = false;
	
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "distortion")
	  distortion = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "yield") {
	  if (utils::ipiece(piter->second) == "source")
	    source = true;
	  else if (utils::ipiece(piter->second) == "target")
	    target = true;
	  else
	    throw std::runtime_error("unknown yield: " + piter->second);
	} else if (utils::ipiece(piter->first) == "goal")
	  goal = piter->second;
	else if (utils::ipiece(piter->first) == "grammar")
	  grammar_local.push_back(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for composer: " << piter->first << "=" << piter->second << std::endl;
      }
	
      if (source && target)
	throw std::runtime_error("Phrase composer can work either source or target yield");
	
      yield_source = source;
    }

    void ComposePhrase::operator()(data_type& data) const
    {
      const lattice_type& lattice = data.lattice;
      hypergraph_type& hypergraph = data.hypergraph;
      hypergraph_type composed;

      hypergraph.clear();
      if (lattice.empty()) return;
    
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      const grammar_type& grammar_compose = (grammar_local.empty() ? grammar : grammar_local);

      utils::resource start;

      grammar_compose.assign(lattice);
    
      cicada::compose_phrase(goal, grammar_compose, distortion, lattice, composed);
    
      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << composed.nodes.size()
		  << " # of edges: " << composed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += composed.nodes.size();
      stat.edge += composed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
    
      hypergraph.swap(composed);
    }

    
    ComposeAlignment::ComposeAlignment(const std::string& parameter,
				       const grammar_type& __grammar,
				       const std::string& __goal,
				       const int __debug)
      : grammar(__grammar),
	goal(__goal),
	lattice_mode(false),
	forest_mode(false),
	debug(__debug)
    { 
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "compose-alignment")
	throw std::runtime_error("this is not a alignment composer");

      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "lattice")
	  lattice_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "forest")
	  forest_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "goal")
	  goal = piter->second;
	else if (utils::ipiece(piter->first) == "grammar")
	  grammar_local.push_back(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for composer: " << piter->first << "=" << piter->second << std::endl;
      }
	
      if (lattice_mode && forest_mode)
	throw std::runtime_error("either lattice or forest");

      if (! lattice_mode && ! forest_mode)
	lattice_mode = true;

      name = std::string("compose-alignment") + (forest_mode ? "-forest" : "-lattice");
    }

    void ComposeAlignment::operator()(data_type& data) const
    {
      const lattice_type& lattice = data.lattice;
      hypergraph_type& hypergraph = data.hypergraph;
      
      if (lattice_mode) {
	if (lattice.empty()) {
	  hypergraph.clear();
	  return;
	}
      } else {
	if (! hypergraph.is_valid())
	  return;
      }
      
      lattice_type target;
      if (! data.targets.empty())
	target = lattice_type(data.targets.front());
	
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      hypergraph_type composed;
      
      const grammar_type& grammar_compose = (grammar_local.empty() ? grammar : grammar_local);
	
      utils::resource start;
      
      if (lattice_mode)
	grammar_compose.assign(lattice, target);
      else
	grammar_compose.assign(hypergraph, target);
	
      if (lattice_mode)
	cicada::compose_alignment(goal, grammar_compose, lattice, target, composed);
      else
	cicada::compose_alignment(goal, grammar_compose, hypergraph, target, composed);
    
      utils::resource end;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << composed.nodes.size()
		  << " # of edges: " << composed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += composed.nodes.size();
      stat.edge += composed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
    
      hypergraph.swap(composed);
    }

    ComposeDependency::ComposeDependency(const std::string& parameter,
					 const int __debug)
      : arc_standard(false),
	arc_eager(false),
	hybrid(false),
	top_down(false),
	degree2(false),
	debug(__debug)
    {
      typedef cicada::Parameter param_type;
	
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "compose-dependency")
	throw std::runtime_error("this is not a dependency composer");
      
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "arc-standard" || utils::ipiece(piter->first) == "standard")
	  arc_standard = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "arc-eager" || utils::ipiece(piter->first) == "eager")
	  arc_eager = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "hybrid")
	  hybrid = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "top-down")
	  top_down = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "degree2")
	  degree2 = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for dependency composer: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (int(arc_standard) + arc_eager + hybrid + top_down + degree2 == 0)
	arc_standard = true;
      
      if (int(arc_standard) + arc_eager + hybrid + top_down + degree2 > 1)
	throw std::runtime_error("you can specify either arc-standard, arc-eager, hybrid, top-down or degree2");
      
      name = std::string("compose-dependency-") + (degree2 ? "degree2" : (hybrid ? "hybrid" : (top_down ? "top-down" : (arc_standard ? "arc-standard" : "arc-eager"))));
    }
    
    void ComposeDependency::operator()(data_type& data) const
    {
      const lattice_type& lattice = data.lattice;
      hypergraph_type& hypergraph = data.hypergraph;
      
      if (lattice.empty()) {
	hypergraph.clear();
	return;
      }
      
      if (debug)
	std::cerr << name << ": " << data.id << std::endl;
      
      hypergraph_type composed;
      
      utils::resource start;
      
      if (hybrid)
	cicada::compose_dependency_hybrid(lattice, composed);
      else if (top_down)
	cicada::compose_dependency_top_down(lattice, composed);
      else if (arc_standard)
	cicada::compose_dependency_arc_standard(lattice, composed);
      else if (arc_eager)
	cicada::compose_dependency_arc_eager(lattice, composed);
      else if (degree2)
	cicada::compose_dependency_degree2(lattice, composed);
      else
	throw std::runtime_error("not implemented?");
      
      utils::resource end;
      
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << std::endl;
    
      if (debug)
	std::cerr << name << ": " << data.id
		  << " # of nodes: " << composed.nodes.size()
		  << " # of edges: " << composed.edges.size()
		  << " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		  << std::endl;

      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.node += composed.nodes.size();
      stat.edge += composed.edges.size();
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
      
      hypergraph.swap(composed);
    }
    
  };
};
