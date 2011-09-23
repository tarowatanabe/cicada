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
    typedef std::vector<coverage_id_type, std::allocator<coverage_id_type> > coverage_id_set_type;
    
    struct Candidate
    {
      typedef boost::array<hypegraph_type::id_type, 2> tails_type;

      score_type score;

      feature_set_type features;
      coverage_id_type coverage;
      tails_type tails;
      
      typename phrase_candidate_set_type::const_iterator phrase_first;
      typename phrase_candidate_set_type::const_iterator phrase_last;
      
      int first;
      int last;
      
      Candidate() : score(), features(), coverage(), tails() {}
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
      
      heaps.clear();
      coverages.clear();
      nodes.clear();
      
      heaps.resize(lattice.size() + 1);
      
      coverage_type __coverage_goal;
      for (size_t i = 0; i != lattice.size(); ++ i)
	__coverage_goal.set(i);
      
      const coverage_id_type coverage_start = coverage_map(coverage_type());
      const coverage_id_type coverage_goal  = coverage_map(__coverage_goal);
      
      for (size_type cardinality = 0; i <= lattice.size(); ++ cardinality) {
	// states...
	// we will synchronize by the candidate_type
	
	candidate_heap_type& heap = heaps[cardinality];
	
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

	  hypergraph_type::id_type node_head = nodes[item->coverage];
	  
	  if (item->tails.front() != hypegraph_type::invalid) {
	    if (item->tails.back() == hypergraph_type::invalid) {
	      item->tails.back() = graph.add_node().id;
	      
	      hypergraph_type::edge_type& edge = grah.add_edge(item->tails.begin(), item->tails.end());
	      edge.rule = rule_x1_x2;
	      
	      graph.connect_edge(edge.id, nodes[item->coverage]);
	    }
	    
	    node_head = item->tails.back();
	  }
	  
	  hypergraph_type::edge_type& edge = graph.add_edge();
	  edge.rule = phrase_cand.rule;
	  edge.features = phrase_cand.features + item->features;
	  edge.attributes = phrase_cand.attributes;
	  
	  edge.attributes[attr_phrase_span_first] = attribute_set_type::int_type(item->first);
	  edge.attributes[attr_phrase_span_last]  = attribute_set_type::int_type(item->last);
	  
	  graph.connect_edge(edge.id, node_head);
	  
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
	coverage_id_set_type::const_iterator citer_end = coverages_cardinality.end();
	for (coverage_id_set_type::const_iterator citer = coverages_cardinality.begin(); citer != citer_end; ++ citer) {
	  typedef std::vector<int, std::allocator<int> > node_set_type;
	  typedef std::pair<transducder_type::id_type, feature_set_type> intersected_type;
	  typedef std::deque<intersected_type, std::allocator<intersected_type> > intersected_set_type;
	  typedef std::vector<intersected_set_type, std::allocator<intersected_set_type> > intersected_map_type;

	  const coverage_type& coverage = coverages[*citer];
	  
	  const int first = coverage.select(1, false);
	  const int last  = utils::bithack::min(static_cast<int>(lattices.size()), first + max_distortion + 1);
	  
	  node_set_type        nodes;
	  node_set_type        nodes_next;
	  coverage_type        visited;
	  intersected_map_type intersected(lattices.size() + 1);
	  
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
		  
		  intersected.clear();
		  intersected.resize(lattice.size() + 1);
		  intersected[*niter].push_back(std::make_pair(transducer.root(), feature_set_type()));
		  
		  const size_type first = *niter;
		  for (size_type last = first; last <= lattice.size() && (first == last || ! coverage.test(last - 1)); ++ last) {
		    if (first != last) {
		      coverage_new.set(last - 1);
		      ++ cardinality_new;
		    }
		    
		    const coverage_id_type coverage_new_id = coverage_map(coverage_new);
		    
		    intersected_set_type::const_iterator niter_end = intersected[last].end();
		    for (intersected_set_type::const_iterator niter = intersected[last].begin(); niter != niter_end; ++ niter) {
		      const phrase_candidate_set_type& phrases = candidate_phrases(table, niter->first);
		      
		      if (! phrases.empty()) {
		      
			candidate_type& cand = candidates.push_back(candidate_type());
			
			cand.score = function(niter->second);
			
			cand.features = niter->second;
			cand.coverage = coverage_new;
			cand.tails.front() = nodes[coverage_new];
			cand.tails.back()  = hypergraph_type::invalid;
			
			cand.phrase_first = phrases.begin();
			cand.phrase_last  = phrasees.end();
			
			cand.first = first;
			cand.last  = last;
			
			heaps[cardinality_new].push(&cand);
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
      
      // goal...
      if (nodes[coverage_goal] != hypergraph_type::invalid) {
	hypergraph_type::edge_type& edge = graph.add_edge(&(nodes[coverage_goal]), &(nodes[coverage_goal]) + 1);
	edge.rule = rule_goal;
	
	edge.attributes[attr_phrase_span_first] = attribute_set_type::int_type(0);
	edge.attributes[attr_phrase_span_last]  = attribute_set_type::int_type(lattice.size());
	
	hypergraph_type::node_type& node = graph.add_node();
	
	graph.connect_edge(edge.id, node.id);
	
	graph.goal = node.id;
	graph.topologically_sort();
      }
    }

    coverage_id_type coverage_map(const coverage_type& coverage)
    {
      coverage_set_type::iterator citer = coverages.insert(coverage).first;
      nodes.resize(coverages.size(), hypergraph_type::invalid);
      return coverages.begin() - citer;
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
