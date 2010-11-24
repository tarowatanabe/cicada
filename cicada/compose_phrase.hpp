// -*- mode: c++ -*-

#ifndef __CICADA__COMPOSE_PHRASE__HPP__
#define __CICADA__COMPOSE_PHRASE__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/bit_vector.hpp>

#include <utils/sgi_hash_set.hpp>
#include <utils/sgi_hash_map.hpp>

namespace cicada
{
  // implementation of 
  //
  // @InProceedings{huang-chiang:2007:ACLMain,
  //  author    = {Huang, Liang  and  Chiang, David},
  //  title     = {Forest Rescoring: Faster Decoding with Integrated Language Models},
  //  booktitle = {Proceedings of the 45th Annual Meeting of the Association of Computational Linguistics},
  //  month     = {June},
  //  year      = {2007},
  //  address   = {Prague, Czech Republic},
  //  publisher = {Association for Computational Linguistics},
  //  pages     = {144--151},
  //  url       = {http://www.aclweb.org/anthology/P07-1019}
  //  }
  //

  struct ComposePhrase
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;

    typedef Lattice    lattice_type;
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;
    
    
    typedef utils::bit_vector<1024> coverage_type;

    struct State
    {
      const coverage_type* coverage;
      
      int grammar_id;
      transducer_type::id_type node;
      
      int first;
      int last;

      feature_set_type features;
      
      State(const coverage_type* __coverage,
	    const int& __grammar_id, const transducer_type::id_type& __node,
	    const int& __first, const int& __last,
	    const feature_set_type& __features)
	: coverage(__coverage),
	  grammar_id(__grammar_id), node(__node),
	  first(__first), last(__last),
	  features(__features) {}
    };

    typedef State state_type;
    
    typedef std::deque<state_type, std::allocator<state_type> > queue_type;

#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<const coverage_type*, hypergraph_type::id_type, boost::hash<const coverage_type*>, std::equal_to<const coverage_type*>,
				    std::allocator<std::pair<const coverage_type*, hypergraph_type::id_type> > > node_map_type;
#else
    typedef sgi::hash_map<const coverage_type*, hypergraph_type::id_type, boost::hash<const coverage_type*>, std::equal_to<const coverage_type*>,
			  std::allocator<std::pair<const coverage_type*, hypergraph_type::id_type> > > node_map_type;
#endif
#ifdef HAVE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_set<coverage_type, boost::hash<coverage_type>, std::equal_to<coverage_type>,
				    std::allocator<coverage_type > > coverage_set_type;
#else
    typedef sgi::hash_set<coverage_type, boost::hash<coverage_type>, std::equal_to<coverage_type>,
			  std::allocator<coverage_type > > coverage_set_type;
#endif
    
    ComposePhrase(const symbol_type& non_terminal,
		  const grammar_type& __grammar,
		  const int& __max_distortion,
		  const bool __yield_source)
      : grammar(__grammar),
	max_distortion(__max_distortion),
	yield_source(__yield_source)
    {
      // initializer...
      rule_goal.reset(new rule_type(vocab_type::GOAL,
				    rule_type::symbol_set_type(1, non_terminal.non_terminal(1))));
      
      std::vector<symbol_type, std::allocator<symbol_type> > sequence(2);
      sequence.front() = non_terminal.non_terminal(1);
      sequence.back()  = non_terminal.non_terminal(2);
      
      rule_x1_x2.reset(new rule_type(non_terminal.non_terminal(),
				     rule_type::symbol_set_type(sequence.begin(), sequence.end())));
    }

    void operator()(const lattice_type& lattice, hypergraph_type& graph)
    {
       graph.clear();

      if (lattice.empty()) return;

      nodes.clear();
      coverages.clear();
      
      queue_type queue;
      
      // breadth first search to construct phrase translational forest
      coverage_type __coverage_goal;
      for (size_t i = 0; i != lattice.size(); ++ i)
	__coverage_goal.set(i);
      
      const coverage_type* coverage_start = coverage_vector(coverage_type()).first;
      const coverage_type* coverage_goal  = coverage_vector(__coverage_goal).first;
      
      nodes[coverage_start] = hypergraph_type::invalid;
      
      // we need to jump the starting positions... and consider the distortion wrt lattice distance..
      // very hard...
      
      const int last = utils::bithack::min(static_cast<int>(lattice.size()), max_distortion + 1);
      
      // breadth first search to compute jump positions...
      {
	typedef std::vector<int, std::allocator<int> > node_set_type;
	
	node_set_type nodes;
	node_set_type nodes_next;
	coverage_type visited;
	
	nodes.push_back(0);
	visited.set(0);
	
	for (int i = 0; i != last && ! nodes.empty(); ++ i) {
	  nodes_next.clear();
	  
	  node_set_type::const_iterator niter_end = nodes.end();
	  for (node_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	    for (size_t table = 0; table != grammar.size(); ++ table)
	      queue.push_back(state_type(coverage_start, table, grammar[table].root(), *niter, *niter, feature_set_type()));
	    
	    lattice_type::arc_set_type::const_iterator aiter_end = lattice[*niter].end();
	    for (lattice_type::arc_set_type::const_iterator aiter = lattice[*niter].begin(); aiter != aiter_end; ++ aiter) {
	      const int next = *niter + aiter->distance;
	      
	      if (! visited[next]) {
		nodes_next.push_back(next);
		visited.set(next);
	      }
	    }
	  }
	  
	  nodes.swap(nodes_next);
	}
      }
      
      while (! queue.empty()) {
	const state_type& state = queue.front();
	
	const transducer_type::rule_pair_set_type& rules = grammar[state.grammar_id].rules(state.node);
	
	if (! rules.empty()) {
	  // next coverage vector...
	  coverage_type __coverage_new = *state.coverage;
	  for (int i = state.first; i != state.last; ++ i)
	    __coverage_new.set(i);

	  std::pair<const coverage_type*, bool> result = coverage_vector(__coverage_new);
	  const coverage_type* coverage_new = result.first;
	  
	  if (result.second) {
	    typedef std::vector<int, std::allocator<int> > node_set_type;
	
	    node_set_type nodes;
	    node_set_type nodes_next;
	    coverage_type visited;
	    
	    const int first = coverage_new->select(1, false);
	    const int last  = utils::bithack::min(static_cast<int>(lattice.size()), first + max_distortion + 1);

	    nodes.push_back(first);
	    visited.set(first);
	    
	    for (int i = first; i != last && ! nodes.empty(); ++ i) {
	      nodes_next.clear();
	      
	      node_set_type::const_iterator niter_end = nodes.end();
	      for (node_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
		if (! coverage_new->test(*niter))
		  for (size_t table = 0; table != grammar.size(); ++ table)
		    queue.push_back(state_type(coverage_new, table, grammar[table].root(), *niter, *niter, feature_set_type()));
		
		lattice_type::arc_set_type::const_iterator aiter_end = lattice[*niter].end();
		for (lattice_type::arc_set_type::const_iterator aiter = lattice[*niter].begin(); aiter != aiter_end; ++ aiter) {
		  const int next = *niter + aiter->distance;
		  
		  if (next < last && ! visited[next]) {
		    nodes_next.push_back(next);
		    visited.set(next);
		  }
		}
	      }

	      nodes.swap(nodes_next);
	    }
	  }
	  
	  // construct graph...
	  
	  hypergraph_type::node_type& node = graph.add_node();

	  transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	  for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	    hypergraph_type::edge_type& edge = graph.add_edge();
	    edge.rule = (yield_source ? riter->source : riter->target);
	    edge.features = riter->features;
	    if (! state.features.empty())
	      edge.features += state.features;
	    
	    edge.first    = state.first;
	    edge.last     = state.last;
	    
	    graph.connect_edge(edge.id, node.id);
	  }
	  
	  node_map_type::const_iterator titer = nodes.find(state.coverage);
	  if (titer == nodes.end())
	    throw std::runtime_error("no precedent nodes?");
	  
	  if (titer->second == hypergraph_type::invalid) {
	    // no precedent nodes... meaning that this is the node created by combining epsilon (or start-coverage)
	    
	    nodes[coverage_new] = node.id;
	  } else {
	    std::pair<node_map_type::iterator, bool> result = nodes.insert(std::make_pair(coverage_new, 0));
	    if (result.second)
	      result.first->second = graph.add_node().id;
	    
	    hypergraph_type::id_type tails[2] = {titer->second, node.id};
	    hypergraph_type::edge_type& edge = graph.add_edge(tails, tails + 2);
	    edge.rule = rule_x1_x2;
	    
	    graph.connect_edge(edge.id, result.first->second);
	  }
	}

	if (state.last != static_cast<int>(lattice.size())) {
	  const lattice_type::arc_set_type& arcs = lattice[state.last];
	  
	  lattice_type::arc_set_type::const_iterator aiter_end = arcs.end();
	  for (lattice_type::arc_set_type::const_iterator aiter = arcs.begin(); aiter != aiter_end; ++ aiter) {
	    const symbol_type& terminal = aiter->label;
	    const int length = aiter->distance;

	    if (terminal == vocab_type::EPSILON) {
	      const size_type rank_first = (state.first == 0 ? size_type(0) : state.coverage->rank(state.first - 1, true));
	      const size_type rank_last  = state.coverage->rank(state.last + length - 1, true);
	      
	      if (rank_first == rank_last)
		queue.push_back(state_type(state.coverage, state.grammar_id, state.node, state.first, state.last + length, state.features + aiter->features));
	      
	    } else {
	      const transducer_type& transducer = grammar[state.grammar_id];
	      const transducer_type::id_type node = transducer.next(state.node, terminal);
	      if (node == transducer.root()) continue;
	      
	      // check if we can use this phrase...
	      // we can check by rank_1(pos) to see whether the number of 1s are equal
	      
	      const size_type rank_first = (state.first == 0 ? size_type(0) : state.coverage->rank(state.first - 1, true));
	      const size_type rank_last  = state.coverage->rank(state.last + length - 1, true);
	      
	      if (rank_first == rank_last)
		queue.push_back(state_type(state.coverage, state.grammar_id, node, state.first, state.last + length, state.features + aiter->features));
	    }
	  }
	}
	
	// finally, pop from the queue!
	queue.pop_front();
      }
      
      node_map_type::iterator niter = nodes.find(coverage_goal);
      if (niter != nodes.end() && niter->second != hypergraph_type::invalid) {
	
	hypergraph_type::edge_type& edge = graph.add_edge(&(niter->second), &(niter->second) + 1);
	edge.rule = rule_goal;
	
	edge.first    = 0;
	edge.last     = lattice.size();
	
	hypergraph_type::node_type& node = graph.add_node();
	
	graph.connect_edge(edge.id, node.id);

	graph.goal = node.id;

	graph.topologically_sort();
      }
    }
    
    std::pair<const coverage_type*, bool> coverage_vector(const coverage_type& coverage)
    {
      std::pair<coverage_set_type::iterator, bool> result = coverages.insert(coverage);
      
      return std::make_pair(&(*result.first), result.second);
    }

  private:
    
    const grammar_type& grammar;
    const int max_distortion;
    const bool yield_source;

    node_map_type     nodes;
    coverage_set_type coverages;
    
    rule_ptr_type rule_goal;
    rule_ptr_type rule_x1_x2;
    
  };
  
  inline
  void compose_phrase(const Symbol& non_terminal, const Grammar& grammar, const Lattice& lattice, const int max_distortion, HyperGraph& graph, const bool yield_source=false)
  {
    ComposePhrase __composer(non_terminal, grammar, max_distortion, yield_source);
    __composer(lattice, graph);
  }

};

#endif
