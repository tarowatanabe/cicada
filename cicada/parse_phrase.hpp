// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_PHRASE__HPP__
#define __CICADA__PARSE_PHRASE__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>


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
  
  template <typename Semiring, typename Function>
  struct ParsePhrase
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;

    typedef Lattice    lattice_type;
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef attribute_set_type::attribute_type attribute_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    typedef Semiring semiring_type;
    typedef Semiring score_type;
    
    typedef Function function_type;
    
    struct PhraseCandidate
    {
      score_type    score;
      
      rule_ptr_type rule;
      
      feature_set_type   features;
      attribute_set_type attributes; 

      PhraseCandidate() : score(), rule(), features(), attributes() {}
      PhraseCandidate() : score(), rule(), features(), attributes() {}
      PhraseCandidate(const score_type& __score, const rule_ptr_type& __rule, const feature_set_type& __features, const attribute_set_type& __attributes)
	: score(__score), rule(__rule), features(__features), attributes(__attributes) {}
      
      void swap(PhraseCandidate& x)
      {
	std::swap(score, x.score);
	rule.swap(x.rule);
	features.swap(x.features);
	attributes.swap(x.attributes);
      }
      
      friend
      void swap(PhraseCandidate& x, PhraseCandidate& y)
      {
	x.swap(y);
      }
    };

    typedef PhraseCandidate phrase_candidate_type;
    typedef utils::simple_vector<phrase_phrase_candidate_type, std::allocator<phrase_candidate_type> > phrase_candidate_set_type;

    
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<transducer_type::id_type, phrase_candidate_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
				    std::allocator<std::pair<const transducer_type::id_type, phrase_candidate_set_type> > > phrase_candidate_map_type;
#else
    typedef sgi::hash_map<transducer_type::id_type, phrase_candidate_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
			  std::allocator<std::pair<const transducer_type::id_type, phrase_candidate_set_type> > > phrase_candidate_map_type;
#endif
    typedef std::vector<phrase_candidate_map_type, std::allocator<phrase_candidate_map_type> > phrase_candidate_table_type;
    
    // coverage vector...
    typedef utils::bit_vector<1024> coverage_type;
    typedef utils::indexed_set<coverage_type, boost::hash<coverage_type>, std::equal_to<coverage_type>, std::allocator<coverage_type> > coverage_set_type;
    typedef coverage_set_type::index_type coverage_id_type;
    
    struct Candidate
    {
      typedef boost::array<hypegraph_type::id_type, 2> tails_type;

      score_type score;
      
      typename phrase_candidate_set_type::const_iterator phrase_first;
      typename phrase_candidate_set_type::const_iterator phrase_last;
      
      feature_set_type features;
      
      coverage_id_type coverage;
      
      tails_type tails;
      
      int first;
      int last;
    };
    
    typedef Candidate candidate_type;
    typedef utils::chunk_vector<candidate_type, 1024 * 8 / sizeof(candidate_type), std::allocator<candidate_type> > candidate_set_type;
    
    struct compare_heap_type
    {
      // we use less, so that when popped from heap, we will grab "greater" in back...
      bool operator()(const candidate_type* x, const candidate_type* y) const
      {
	return x->score * x->first->score < y->score * y->first->score;
      }
    };
    
    typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_base_type;
    typedef utils::std_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type> candidate_heap_type;
    typedef std::vector<candidate_heap_type, std::allocator<candidate_heap_type> > candidate_heap_map_type;
    
    ParsePhrase(const symbol_type& non_terminal,
		const grammar_type& __grammar,
		const function_type& __function,
		const int __beam_size,
		const int& __max_distortion,
		const bool __yield_source)
      : grammar(__grammar),
	function(__function),
	beam_size(__beam_size),
	max_distortion(__max_distortion),
	yield_source(__yield_source),
	attr_phrase_span_first("phrase-span-first"),
	attr_phrase_span_last("phrase-span-last")
    {
      // initializer...
      rule_goal = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, non_terminal.non_terminal(1))));
      
      std::vector<symbol_type, std::allocator<symbol_type> > sequence(2);
      sequence.front() = non_terminal.non_terminal(1);
      sequence.back()  = non_terminal.non_terminal(2);
      
      rule_x1_x2 = rule_type::create(rule_type(non_terminal.non_terminal(), sequence.begin(), sequence.end()));
    }

    void operator()(const lattice_type& lattice, hypergraph_type& graph)
    {
      graph.clear();
      
      if (lattice.empty()) return;
      
      // initialization...
      
      coverages.clear();
      nodes.clear();
      
      
      for (size_type i = 0; i <= lattice.size(); ++ i) {
	// states...
	// we will synchronize by the candidate_type

	const size_type cardinality = i;
	
	candidate_heap_type& heap = heaps[i];
	
	coverages_cardinality.clear();
	if (cardinality == 0)
	  coverages_cardinality.push_back(coverage_start);
	
	for (int num_pop = 0; ! heap.empty() && num_pop != beam_size; /**/) {
	  candidate_type* item = heap.top();
	  heap.pop();

	  const phrase_candidate_type& phrase_cand = *(item->phrase_first);
	  
	  // when constructig hypergraph, we will construct by
	  //
	  // head -> tail tail-for-phrase
	  //
	  // thus, for each candidate, we need to keep head associated with coverage vector + tail-for-phrase associated
	  // with each candidate_type...
	  
	  if (nodes[item->coverage] == hypergraph_type::invalid) {
	    nodes[item->coverage] = graph.add_node().id;
	    coverages_cardinality.push_back(item->coverage);
	  }

	  const hypergraph_type::id_type node_head = nodes[item->coverage];
	  
	  if (item->tails[0] == hypergraph_type::invalid) {
	    // this is the initial phrases...
	    // we will construct a phral edge with head associated with "coverage" node-id
	    
	    // add hyperedge connecting phrase with node_head
	    
	    hypergraph_type::edge_type& edge = graph.add_edge();
	    edge.rule = phrase_cand.rule;
	    edge.features = phrase_cand.features + item->features;
	    edge.attributes = phrase_cand.attributes;
	    
	    edge.attributes[attr_phrase_span_first] = attribute_set_type::int_type(item->first);
	    edge.attributes[attr_phrase_span_last]  = attribute_set_type::int_type(item->last);
	    
	    graph.connect_edge(edge.id, node_head);
	  } else {
	    if (item->tails[1] == hypergraph_type::invalid) {
	      item->tails[1] = graph.add_node().id;
	      
	      // add hyperedge connecting node_head as head, and item->tail and item->tail_phrase as tails.
	      
	      hypergraph_type::edge_type& edge = graph.add_edge();
	      edge.rule = phrase_cand.rule;
	      edge.features = phrase_cand.features + item->features;
	      
	      edge.attributes = phrase_cand.attributes;
	      edge.attributes[attr_phrase_span_first] = attribute_set_type::int_type(item->first);
	      edge.attributes[attr_phrase_span_last]  = attribute_set_type::int_type(item->last);
	      
	      graph.connect_edge(edge.id, item->tail_phrase);
	    }
	    
	    // add hyperedge connecting phrase with item->tail_phrase;
	    hypergraph_type::edge_type& edge = grah.add_edge(item->tails.begin(), item->tails.end());
	    edge.rule = rule_x1_x2;
	    
	    graph.connect_edge(edge.id, node_head);
	  }
	  
	  // proceed to the next...
	  ++ item->phrase_first;
	  if (item->phrase_first != item->phrase_last)
	    heap.push(item);
	  
	  ++ num_pop;
	}
	
	// clear unused heap...
	heap.clear();
	candidate_heap_type(heap).swap(heap);
	
	// enumerate coverages and add new candidates to the heaps...
	coverage_ptr_set_type::const_iterator citer_end = coverages_cardinality.end();
	for (coverage_ptr_set_type::const_iterator citer = coverages_cardinality.begin(); citer != citer_end; ++ citer) {
	  const coverage_type& coverage = coverages[*citer];
	  
	  const int first = coverage.select(1, false);
	  const int last  = utils::bithack::min(static_cast<int>(lattices.size()), first + max_distortion + 1);
	  
	  node_set_type nodes;
	  node_set_type nodes_next;
	  coverage_type visited;
	  
	  nodes.push_back(first);
	  visited.set(first);
	  
	  for (int i = first; i != last && ! nodes.empty(); ++ i) {
	    nodes_next.clear();
	    
	    node_set_type::const_iterator niter_end = nodes.end();
	    for (node_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	      if (! coverage.test(*niter))
		for (size_type table = 0; table != grammar.size(); ++ table) {
		  // start enumerating grammar...
		  
		  const transducer_type& transducer = grammar[table];
		  
		  coverage_type coverage_new    = coverage;
		  size_type     cardinality_new = cardinality;
		  
		  nodes_intersected.clear();
		  nodes_intersected.resize(lattice.size() + 1);
		  nodes_intersected[*niter].push_back(std::make_pair(transducer.root(), feature_set_type()));
		  
		  const size_type first = *niter;
		  for (size_type last = first; last <= lattice.size() && (first == last || ! coverage.test(last - 1)); ++ last) {
		    if (first != last) {
		      coverage_new.set(last - 1);
		      ++ cardinality_new;
		    }
		    
		    node_intersected_set_type::const_iterator niter_end = nodes_intersected[last].end();
		    for (node_intersected_set_type::const_iterator niter = nodes_intersected[last].begin(); niter != niter_end; ++ niter) {
		      const phrase_candidate_set_type& phrases = candidate_phrases(table, niter->first);
		      
		      if (phrases.empty()) continue;
		      
		      candidate_type& cand = candidates.push_back(candidate_type());
		      
		      cand.score = function(niter->second);
		      
		      cand.phrase_first = phrases.begin();
		      cand.phrase_last  = phrasees.end();
		      
		      cand.features = niter->second;
		      
		      cand.tail = ; // this is the node-id associated with the previous coverage...
		      cand.tail_phrase = hypergraph_type::invalid;
		      
		      cand.first = first;
		      cand.last  = last;
		      
		      heaps[cardinality_new].push(&cand);
		      
		      // extention...
		      if (last != lattice.size()) {
			const lattice_type::arc_set_type& arcs = lattice[last];
			
			lattice_type::arc_set_type::const_iterator aiter_end = arcs.end();
			for (lattice_type::arc_set_type::const_iterator aiter = arcs.begin(); aiter != aiter_end; ++ aiter) {
			  const symbol_type& terminal = aiter->label;
			  const int length = aiter->distance;
			  
			  if (terminal == vocab_type::EPSILON)
			    nodes_intersected[last + length].push_back(std::make_pair(niter->first, niter->second + aiter->features));
			  else {
			    const transducer_type::id_type node = transducer.next(niter->first, terminal);
			    if (node == transducer.root()) continue;
			    
			    nodes_intersected[last + length].push_back(std::make_pair(node, niter->second + aiter->features));
			  }
			}
		      }
		    }
		  }
		}
	      
	      // fill-in nodes-next to determine the jump position over the lattice...
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
      }
      
      // goal...
      
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
	      for (node_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter)
		if (! coverage_new->test(*niter)) {
		  for (size_t table = 0; table != grammar.size(); ++ table)
		    queue.push_back(state_type(coverage_new, table, grammar[table].root(), *niter, *niter, feature_set_type()));
		  
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
	  
	  // construct graph...
	  
	  hypergraph_type::node_type& node = graph.add_node();

	  transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	  for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	    hypergraph_type::edge_type& edge = graph.add_edge();
	    edge.rule = (yield_source ? riter->source : riter->target);
	    edge.features = riter->features;
	    if (! state.features.empty())
	      edge.features += state.features;
	    edge.attributes = riter->attributes;
	    
	    // assign metadata...
	    edge.attributes[attr_phrase_span_first] = attribute_set_type::int_type(state.first);
	    edge.attributes[attr_phrase_span_last]  = attribute_set_type::int_type(state.last);
	    
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
	
	edge.attributes[attr_phrase_span_first] = attribute_set_type::int_type(0);
	edge.attributes[attr_phrase_span_last]  = attribute_set_type::int_type(lattice.size());
	
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
    
    
    template <typename Tp>
    struct greater_score
    {
      bool operator()(const Tp& x, const Tp& y) const
      {
	return x.score > y.score;
      }
    };
    
    const phrase_candidate_set_type& candidate_phrases(const size_type& table, const transducer_type::id_type& node)
    {
      typename phrase_candidate_map_type::iterator riter = phrase_tables[table].find(node);
      if (riter == phrase_tables[table].end()) {
	const transducer_type::rule_pair_set_type& phrases = grammar[table].rules(node);
	
	riter = phrase_tables[table].insert(std::make_pair(node, phrase_candidate_set_type(phrases.size()))).first;
	
	typename phrase_candidate_set_type::iterator citer = riter->second.begin();
	transducer_type::rule_pair_set_type::const_iterator iter_end = phrases.end();
	for (transducer_type::rule_pair_set_type::const_iterator iter = phrases.begin(); iter != iter_end; ++ iter, ++ citer)
	  *citer = phrase_candidate_type(function(iter->features),
					 yield_source ? iter->source : iter->target,
					 iter->features,
					 iter->attributes);
	
	std::sort(riter->second.begin(), riter->second.end(), greater_score<phrase_candidate_type>());
      }
      
      return riter->second;
    }

  private:
    
    const grammar_type& grammar;
    const function_type& function;
    const int beam_size;
    
    const int max_distortion;
    const bool yield_source;
    
    const attribute_type attr_phrase_span_first;
    const attribute_type attr_phrase_span_last;

    node_map_type     nodes;
    coverage_set_type coverages;
    
    rule_ptr_type rule_goal;
    rule_ptr_type rule_x1_x2;
  };
  
  template <typename Function>
  inline
  void parse_phrase(const Symbol& non_terminal, const Grammar& grammar, const Function& function, const int beam_size, const int max_distortion, const Lattice& lattice, HyperGraph& graph, const bool yield_source=false)
  {
    ParsePhrase<typename FUnction::value_type, Function> __parser(non_terminal, grammar, function, beam_size, max_distortion, yield_source);
    __parser(lattice, graph);
  }

};

#endif
