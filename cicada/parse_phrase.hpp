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
    typedef utils::simple_vector<phrase_candidate_type, std::allocator<phrase_candidate_type> > phrase_candidate_set_type;

    
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
    typedef std::vector<coverage_id_type, std::allocator<coverage_id_type> > coverage_id_set_type;
    
    struct LatticeEdge
    {
      LatticeEdge(const coverage_id_type& __tail, const coverage_id_type& __head)
	: tail(__tail), head(__head) {}
      LatticeEdge() {}
      
      coverage_id_type tail;
      coverage_id_type head;

      score_type       score;
      feature_set_type features;
      
      typename phrase_candidate_set_type::const_iterator phrase_first;
      typename phrase_candidate_set_type::const_iterator phrase_last;
      
      int first;
      int last;
    };
    typedef LatticeEdge lattice_edge_type;
    typedef utils::chunk_vector<lattice_edge_type, 4096 / sizeof(lattice_edge_type), std::allocator<lattice_edge_type> > lattice_edge_set_type;
    
    struct LatticeNode
    {
      LatticeNode() {}
      LatticeNode(const coverage_id_type& __coverage) : coverage(__coverage) {}

      coverage_id_type      coverage;
      lattice_edge_set_type edges;
    };
    typedef LatticeNode lattice_node_type;
    typedef utils::chunk_vector<lattice_node_type, 4096 / sizeof(lattice_node_type), std::allocator<lattice_node_type> > lattice_node_set_type;
    typedef std::vector<lattice_node_set_type, std::allocator<lattice_node_set_type> > lattice_node_map_type;
    
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_map_type;
    typedef std::vector<score_type,  std::allocator<score_type> >                            score_map_type;
    
    struct Candidate
    {
      const lattice_edge_type* edge;
      
      typename phrase_candidate_set_type::const_iterator first;
      typename phrase_candidate_set_type::const_iterator last;
      
      score_type score;
      hypergraph_type::id_type tail;
      
      Candidate() : edge(0) {}
      Candidate(const lattice_edge_type& __edge)
	: edge(&__edge),
	  first(__edge.phrase_first),
	  last(__edge.phrase_last),
	  score(__edge.score),
	  tail(hypergraph_type::invalid) {}
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
      lattice_nodes.clear();
      coverages.clear();
      
      node_map_forward.clear();
      node_map_backward.clear();
      forward.clear();
      backward.clear();
      
      candidates.clear();
      heaps.clear();

      phrase_tables.clear();
      phrase_tables.resize(grammar.size());
      
      const coverage_id_type coverage_start_id = coverage_map(coverage_type());
      
      // forward computation: we will construct a "lazy" lattice
      
      lattice_nodes.resize(lattice.size() + 1);
      forward[coverage_start_id] = semiring::traits<score_type>::one();
      node_map_forward[coverage_start_id] = lattice_nodes[0].size();
      lattice_nodes[0].push_back(lattice_node_type(coverage_start_id));
      
      for (size_type cardinality = 0; cardinality <= lattice.size(); ++ cardinality) {
	
	// enumerate coverages and add new candidates to the heaps...
	typename lattice_node_set_type::const_iterator citer_end = lattice_nodes[cardinality].end();
	for (typename lattice_node_set_type::const_iterator citer = lattice_nodes[cardinality].begin(); citer != citer_end; ++ citer) {
	  typedef std::vector<int, std::allocator<int> > node_set_type;
	  typedef std::pair<transducer_type::id_type, feature_set_type> intersected_type;
	  typedef std::deque<intersected_type, std::allocator<intersected_type> > intersected_set_type;
	  typedef std::vector<intersected_set_type, std::allocator<intersected_set_type> > intersected_map_type;
	  
	  const coverage_id_type& coverage_id = citer->coverage;
	  
	  const int first = coverages[coverage_id].select(1, false);
	  const int last  = utils::bithack::min(static_cast<int>(lattice.size()), first + max_distortion + 1);
	  
	  node_set_type        nodes;
	  node_set_type        nodes_next;
	  coverage_type        visited;
	  intersected_map_type intersected(lattice.size() + 1);

	  nodes.push_back(first);
	  visited.set(first);
	  
	  for (int i = first; i != last && ! nodes.empty(); ++ i) {
	    nodes_next.clear();

	    node_set_type::const_iterator niter_end = nodes.end();
	    for (node_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	      if (! coverages[coverage_id].test(*niter))
		for (size_type table = 0; table != grammar.size(); ++ table) {
		  // start enumerating grammar...
		  
		  const transducer_type& transducer = grammar[table];
		  
		  coverage_type coverage_new    = coverages[coverage_id];
		  size_type     cardinality_new = cardinality;
		  
		  intersected.clear();
		  intersected.resize(lattice.size() + 1);
		  intersected[*niter].push_back(std::make_pair(transducer.root(), feature_set_type()));
		  
		  const size_type first = *niter;
		  for (size_type last = first; last <= lattice.size() && (first == last || ! coverages[coverage_id].test(last - 1)); ++ last) {
		    if (first != last) {
		      coverage_new.set(last - 1);
		      ++ cardinality_new;
		    }
		    
		    const coverage_id_type coverage_new_id = coverage_map(coverage_new);
		    
		    if (node_map_forward[coverage_new_id] == hypergraph_type::invalid) {
		      node_map_forward[coverage_new_id] = lattice_nodes[cardinality_new].size();
		      lattice_nodes[cardinality_new].push_back(lattice_node_type(coverage_new_id));
		    }

		    lattice_edge_set_type& lattice_edges = lattice_nodes[cardinality_new][node_map_forward[coverage_new_id]].edges;
		    
		    intersected_set_type::const_iterator niter_end = intersected[last].end();
		    for (intersected_set_type::const_iterator niter = intersected[last].begin(); niter != niter_end; ++ niter) {
		      const phrase_candidate_set_type& phrases = candidate_phrases(table, niter->first);
		      
		      if (! phrases.empty()) {
			lattice_edges.push_back(lattice_edge_type(coverage_id, coverage_new_id));
			lattice_edge_type& lattice_edge = lattice_edges.back();
			
			lattice_edge.score    = function(niter->second);
			lattice_edge.features = niter->second;
			
			lattice_edge.phrase_first = phrases.begin();
			lattice_edge.phrase_last  = phrases.end();
			
			lattice_edge.first = first;
			lattice_edge.last  = last;
			
			const score_type score = forward[coverage_id] * lattice_edge.score * lattice_edge.phrase_first->score;
			
			forward[coverage_new_id] = std::max(forward[coverage_new_id], score);
		      }
		      
		      // extention...
		      if (last != lattice.size()) {
			const lattice_type::arc_set_type& arcs = lattice[last];
			
			lattice_type::arc_set_type::const_iterator aiter_end = arcs.end();
			for (lattice_type::arc_set_type::const_iterator aiter = arcs.begin(); aiter != aiter_end; ++ aiter) {
			  const symbol_type& terminal = aiter->label;
			  const int length = aiter->distance;
			  
			  if (terminal == vocab_type::EPSILON)
			    intersected[last + length].push_back(std::make_pair(niter->first, niter->second + aiter->features));
			  else {
			    const transducer_type::id_type node = transducer.next(niter->first, terminal);
			    if (node == transducer.root()) continue;

			    intersected[last + length].push_back(std::make_pair(node, niter->second + aiter->features));
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
      
      // check whether we have reached a goal...
      if (lattice_nodes.back().empty()) return;
      if (lattice_nodes.back().size() != 1)
	throw std::runtime_error("invalid goal node...?");
      
      // we will start backward path, which enumerate edges with "beam"
      
      const coverage_id_type coverage_goal_id = lattice_nodes.back().front().coverage;
      
      heaps.resize(lattice.size() + 1);
      node_map_backward.resize(coverages.size(), hypergraph_type::invalid);
      backward.resize(coverages.size());
      
      backward[coverage_goal_id] = semiring::traits<score_type>::one();
      node_map_backward[coverage_goal_id] = graph.add_node().id;
      
      coverage_id_set_type coverages_backward;
      coverages_backward.push_back(coverage_goal_id);
      
      for (int cardinality = lattice.size(); cardinality >= 0; -- cardinality) {
	// enumerate heaps...

	// when cardinality == lattice.size(), we will simply skip this, since nothing is pushed into headp
	candidate_heap_type& heap = heaps[cardinality];
	for (int num_pop = 0; ! heap.empty() && num_pop != beam_size; /**/) {
	  const candidate_type* item = heap.top();
	  heap.pop();
	  
	  const lattice_edge_type& lattice_edge = *(item->edge);
	  
	  // update backward score...
	  backward[lattice_edge.tail] = std::max(backward[lattice_edge.tail], item->score * item->first->score / forward[lattice_edge.tail]);
	  
	  hypergraph_type::id_type head = node_map_backward[lattice_edge.head];

	  if (head == hypergraph_type::invalid)
	    throw std::runtime_error("invalid head node id");
	  
	  if (lattice_edge.tail != coverage_start_id) {
	    if (node_map_backward[lattice_edge.tail] == hypergraph_type::invalid) {
	      node_map_backward[lattice_edge.tail] = graph.add_node().id;
	      coverages_backward.push_back(lattice_edge.tail);
	    }
	    
	    if (item->tail == hypergraph_type::invalid) {
	      const_cast<candidate_type*>(item)->tail = graph.add_node().id;
	      
	      const hypergraph_type::id_type tails[2] = {node_map_backward[lattice_edge.tail], item->tail};
	      
	      hypergraph_type::edge_type& edge = graph.add_edge(tails, tails + 2);
	      edge.rule = rule_x1_x2;
	      
	      graph.connect_edge(edge.id, node_map_backward[lattice_edge.head]);
	    }
	    
	    head = item->tail;
	  }
	  
	  hypergraph_type::edge_type& edge = graph.add_edge();
	  edge.rule = item->first->rule;
	  edge.features = lattice_edge.features + item->first->features;
	  edge.attributes = item->first->attributes;
	  
	  edge.attributes[attr_phrase_span_first] = attribute_set_type::int_type(lattice_edge.first);
	  edge.attributes[attr_phrase_span_last]  = attribute_set_type::int_type(lattice_edge.last);
	  
	  graph.connect_edge(edge.id, head);

	  ++ const_cast<candidate_type*>(item)->first;
	  if (item->first != item->last)
	    heap.push(item);
	}
	heap.clear();
	
	// enumerate next coverages.
	// coverages_backward is updated during initialization or the items in the heaps.
	coverage_id_set_type::const_iterator citer_end = coverages_backward.end();
	for (coverage_id_set_type::const_iterator citer = coverages_backward.begin(); citer != citer_end; ++ citer) {
	  const coverage_id_type& coverage_id = *citer;
	  const coverage_type& coverage = coverages[coverage_id];
	  const int cardinality = coverage.count();
	  
	  const lattice_edge_set_type& lattice_edges = lattice_nodes[cardinality][node_map_forward[coverage_id]].edges;
	  
	  typename lattice_edge_set_type::const_iterator eiter_end = lattice_edges.end();
	  for (typename lattice_edge_set_type::const_iterator eiter = lattice_edges.begin(); eiter != eiter_end; ++ eiter) {
	    const lattice_edge_type& lattice_edge = *eiter;
	    const int cardinality_prev = cardinality - (lattice_edge.last - lattice_edge.first);

	    candidates.push_back(candidate_type(lattice_edge));
	    candidate_type& cand = candidates.back();
	    
	    // for comparision, we need to integrate forward estimates
	    cand.score = forward[lattice_edge.tail] * lattice_edge.score * backward[lattice_edge.head];
	    
	    heaps[cardinality_prev].push(&cand);
	  }
	}
	
	// finished.. chear here...
	coverages_backward.clear();
      }
      
      // add extra rule_goal hyperedge
      hypergraph_type::edge_type& edge = graph.add_edge(&node_map_backward[coverage_goal_id], (&node_map_backward[coverage_goal_id]) + 1);
      edge.rule = rule_goal;
      
      graph.goal = graph.add_node().id;
      graph.connect_edge(edge.id, graph.goal);
      
      graph.topologically_sort();
    }
    
    coverage_id_type coverage_map(const coverage_type& coverage)
    {
      coverage_set_type::iterator citer = coverages.insert(coverage).first;
      node_map_forward.resize(coverages.size(), hypergraph_type::invalid);
      forward.resize(coverages.size(), score_type());
      return citer - coverages.begin();
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

    coverage_set_type     coverages;
    lattice_node_map_type lattice_nodes;
    
    node_map_type node_map_forward;
    node_map_type node_map_backward;
    
    score_map_type forward;
    score_map_type backward;
    
    candidate_set_type      candidates;
    candidate_heap_map_type heaps;

    phrase_candidate_table_type phrase_tables;
    
    rule_ptr_type rule_goal;
    rule_ptr_type rule_x1_x2;
  };
  
  template <typename Function>
  inline
  void parse_phrase(const Symbol& non_terminal, const Grammar& grammar, const Function& function, const int beam_size, const int max_distortion, const Lattice& lattice, HyperGraph& graph, const bool yield_source=false)
  {
    ParsePhrase<typename Function::value_type, Function> __parser(non_terminal, grammar, function, beam_size, max_distortion, yield_source);
    __parser(lattice, graph);
  }

};

#endif
