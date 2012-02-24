// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_AGENDA__HPP__
#define __CICADA__PARSE_AGENDA__HPP__ 1

#include <vector>
#include <algorithm>
#include <set>
#include <iostream>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/b_heap.hpp>
#include <utils/std_heap.hpp>
#include <utils/bithack.hpp>
#include <utils/dense_hash_map.hpp>
#include <utils/dense_hash_set.hpp>

namespace cicada
{
  
  // an agenda-based chart parsing algorithm processing in bottom-up fashion
  //
  // @INPROCEEDINGS{Klein01parsingand,
  //  author = {Dan Klein and Christopher D. Manning},
  //  title = {Parsing and Hypergraphs},
  //  booktitle = {IN IWPT},
  //  year = {2001},
  //  pages = {123--134},
  //  publisher = {}
  // }
  // 
  // One of the main reason for "bottom-up" is the grammar encoding represented in transducer class
  // which do not support Earlye-style lhs-first encoding and differentiating terminals/non-terminals
  //

  template <typename Semiring, typename Function>
  struct ParseAgenda
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

    struct RuleCandidate
    {
      score_type    score;
      
      rule_ptr_type rule;
      
      feature_set_type   features;
      attribute_set_type attributes;
      
      RuleCandidate() : score(), rule(), features(), attributes() {}
      RuleCandidate(const score_type& __score, const rule_ptr_type& __rule, const feature_set_type& __features, const attribute_set_type& __attributes)
	: score(__score), rule(__rule), features(__features), attributes(__attributes) {}
    };
    typedef RuleCandidate rule_candidate_type;
    typedef utils::chunk_vector<rule_candidate_type, 4096 / sizeof(rule_candidate_type), std::allocator<rule_candidate_type> > rule_candidate_set_type;
    
    
    typedef std::vector<const rule_candidate_type*, std::allocator<const rule_candidate_type*> > rule_candidate_ptr_set_type;
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<transducer_type::id_type, rule_candidate_ptr_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
				    std::allocator<std::pair<const transducer_type::id_type, rule_candidate_ptr_set_type> > > rule_candidate_map_type;
#else
    typedef sgi::hash_map<transducer_type::id_type, rule_candidate_ptr_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
			  std::allocator<std::pair<const transducer_type::id_type, rule_candidate_ptr_set_type> > > rule_candidate_map_type;
#endif
    typedef std::vector<rule_candidate_map_type, std::allocator<rule_candidate_map_type> > rule_candidate_table_type;
    
    // since we need to differentiate by lhs, we need to check "pos"
    struct Dot
    {
      size_type                table;
      transducer_type::id_type node;
      
      Dot() : table(size_type(-1)), node() {}
      Dot(const size_type& __table, const transducer_type::id_type& __node)
	: table(__table), node(__node) {}

      friend
      bool operator==(const Dot& x, const Dot& y)
      {
	return x.table == y.table && x.node == y.node;
      }
      
      friend
      std::ostream& operator<<(std::ostream& os, const Dot& x)
      {
	os << x.table << ".." << x.node;
	return os;
      }
      
    };
    
    typedef Dot dot_type;

    struct Span
    {
      int first;
      int last;
      int level;
      
      Span() : first(0), last(0), level(0) {}
      Span(const int& __first, const int& __last) : first(__first), last(__last), level(0) {}
      Span(const int& __first, const int& __last, const int& __level) : first(__first), last(__last), level(__level) {}

      size_type size() const { return last - first; }

      friend
      bool operator==(const Span& x, const Span& y)
      {
	return x.first == y.first && x.last == y.last && x.level == y.level;
      }
      
      friend
      std::ostream& operator<<(std::ostream& os, const Span& x)
      {
	os << x.first << ".." << x.last << ':' << x.level;
	return os;
      }
    };

    typedef Span span_type;

    struct Edge
    {
      typedef Edge edge_type;
      
      score_type  score;
      
      const edge_type* active;
      const edge_type* passive;
      
      dot_type    dot;
      span_type   span;
      
      const rule_candidate_type* rule;
      feature_set_type           features;
      attribute_set_type         attributes;

      // introduction edge
      Edge(const score_type& __score,
	   const dot_type& __dot,
	   const span_type& __span)
	: score(__score),
	  active(0), passive(0), 
	  dot(__dot), span(__span),
	  rule(0) {}
      
      // scanning
      Edge(const score_type& __score,
	   const edge_type& __active,
	   const dot_type& __dot,
	   const span_type& __span,
	   const feature_set_type& __features, 
	   const attribute_set_type& __attributes)
	: score(__score),
	  active(&__active), passive(0), 
	  dot(__dot), span(__span),
	  rule(0),
	  features(__features),
	  attributes(__attributes) {}
      Edge(const score_type& __score,
	   const edge_type& __active,
	   const dot_type& __dot,
	   const span_type& __span,
	   const rule_candidate_type* __rule,
	   const feature_set_type& __features, 
	   const attribute_set_type& __attributes)
	: score(__score),
	  active(&__active), passive(0), 
	  dot(__dot), span(__span),
	  rule(__rule),
	  features(__features),
	  attributes(__attributes) {}
      
      // prediction
      Edge(const score_type& __score,
	   const edge_type& __passive,
	   const dot_type& __dot,
	   const span_type& __span)
	: score(__score),
	  active(0), passive(&__passive), 
	  dot(__dot), span(__span),
	  rule(0) {}
      Edge(const score_type& __score,
	   const edge_type& __passive,
	   const dot_type& __dot,
	   const span_type& __span,
	   const rule_candidate_type* __rule)
	: score(__score),
	  active(0), passive(&__passive), 
	  dot(__dot), span(__span),
	  rule(__rule) {}
      
      // completion
      Edge(const score_type& __score,
	   const edge_type& __active,
	   const edge_type& __passive,
	   const dot_type& __dot,
	   const span_type& __span,
	   const feature_set_type& __features, 
	   const attribute_set_type& __attributes)
	: score(__score),
	  active(&__active), passive(&__passive), 
	  dot(__dot), span(__span),
	  rule(0),
	  features(__features),
	  attributes(__attributes) {}
      Edge(const score_type& __score,
	   const edge_type& __active,
	   const edge_type& __passive,
	   const dot_type& __dot,
	   const span_type& __span,
	   const rule_candidate_type* __rule,
	   const feature_set_type& __features, 
	   const attribute_set_type& __attributes)
	: score(__score),
	  active(&__active), passive(&__passive), 
	  dot(__dot), span(__span),
	  rule(__rule),
	  features(__features),
	  attributes(__attributes) {}
      
    public:
      bool is_passive() const { return rule; }
      bool is_active() const { return ! rule; }
    };
    
    typedef Edge edge_type;
    typedef utils::chunk_vector<edge_type, 1024 * 16 / sizeof(edge_type), std::allocator<edge_type> > edge_set_type;
    typedef std::vector<const edge_type*, std::allocator<const edge_type*> > agenda_exploration_type;
    
    struct edge_heap_type
    {
      bool operator()(const edge_type* x, const edge_type* y) const
      {
	return x->score < y->score;
      }
    };
    
    typedef std::vector<const edge_type*, std::allocator<const edge_type*> > agenda_finishing_base_type;
    typedef utils::std_heap<const edge_type, agenda_finishing_base_type, edge_heap_type> agenda_finishing_type;
    typedef std::vector<agenda_finishing_type, std::allocator<agenda_finishing_type> > agenda_finishing_set_type;
    
    struct Traversal
    {
      const edge_type* active;
      const edge_type* passive;
      bool is_active;
      
      Traversal(const edge_type* __active, const edge_type* __passive, const bool __is_active)
	: active(__active), passive(__passive), is_active(__is_active) {}
      Traversal()
	: active(0), passive(0), is_active(0) {}
    };
    typedef Traversal traversal_type;
    
    typedef utils::hashmurmur<size_type> traversal_hash_type;
    struct traversal_equal_type
    {
      bool operator()(const traversal_type& x, const traversal_type& y) const
      {
	return x.passive == y.passive && x.active == y.active && x.is_active == y.is_active;
      }
    };
    
    typedef google::dense_hash_set<traversal_type, traversal_hash_type, traversal_equal_type > traversal_set_type;
    
    typedef std::vector<const edge_type*, std::allocator<const edge_type*> >   edge_ptr_set_type;
    typedef std::vector<edge_ptr_set_type, std::allocator<edge_ptr_set_type> > edge_set_active_type;
    typedef std::vector<edge_ptr_set_type, std::allocator<edge_ptr_set_type> > edge_set_passive_type;
    
    struct edge_active_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const edge_type* x) const
      {
	return (x ? hasher_type::operator()(x->dot, hasher_type::operator()(x->span)) : size_t(0));
      }
    };

    struct edge_active_equal_type
    {
      bool operator()(const edge_type* x, const edge_type* y) const
      {
	return ((x == y) || (x && y && x->span == y->span && x->dot == y->dot));
      }
    };
    
    struct edge_passive_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const edge_type* x) const
      {
	return (x ? hasher_type::operator()(x->span, x->rule->rule->lhs.id()) : size_t(0));
      }
    };
    struct edge_passive_equal_type
    {
      bool operator()(const edge_type* x, const edge_type* y) const
      {
	return ((x == y) || (x && y && x->span == y->span && x->rule->rule->lhs == y->rule->rule->lhs));
      }
    };
    
    struct head_edge_set_type
    {
      score_type               score;
      hypergraph_type::id_type head;
      edge_ptr_set_type        edges;
      
      head_edge_set_type() : score(), head(hypergraph_type::invalid), edges() {}
    };
    
    typedef google::dense_hash_set<const edge_type*, edge_active_hash_type, edge_active_equal_type > discovered_active_type;
    typedef google::dense_hash_map<const edge_type*, head_edge_set_type, edge_passive_hash_type, edge_passive_equal_type > discovered_passive_type;
    
    
    ParseAgenda(const symbol_type& __goal,
		const grammar_type& __grammar,
		const function_type& __function,
		const int __beam_size,
		const bool __yield_source=false,
		const bool __treebank=false,
		const bool __pos_mode=false)
      : goal(__goal),
	grammar(__grammar),
	function(__function),
	beam_size(__beam_size),
	yield_source(__yield_source),
	treebank(__treebank),
	pos_mode(__pos_mode),
	attr_span_first("span-first"),
	attr_span_last("span-last")
    {
      traversals.set_empty_key(traversal_type());
      discovered_active.set_empty_key(0);
      discovered_passive.set_empty_key(0);
    }

    struct __extract
    {
      const symbol_type& operator()(const symbol_type& x) const
      {
	return x;
      }
    };

    struct __extract_terminal
    {
      symbol_type operator()(const symbol_type& x) const
      {
	return x.terminal();
      }
    };
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      if (lattice.empty()) return;
      
      edges.clear();
      
      agenda_exploration.clear();
      agenda_finishing.clear();
      agenda_finishing.reserve(lattice.size() + 1);
      agenda_finishing.resize(lattice.size() + 1);
      
      traversals.clear();
      
      discovered_active.clear();
      discovered_passive.clear();
      
      edges_active.clear();
      edges_passive.clear();

      edges_active.reserve(lattice.size() + 1);
      edges_passive.reserve(lattice.size() + 1);

      edges_active.resize(lattice.size() + 1);
      edges_passive.resize(lattice.size() + 1);
            
      rule_candidates.clear();
      rule_tables.clear();
      rule_tables.reserve(grammar.size());
      rule_tables.resize(grammar.size());
      
      // initialize by empty word... we do not support epsilon insertion as in the "Parsing and Hypergraphs"
      for (size_type table = 0; table != grammar.size(); ++ table) {
	const transducer_type& transducer = grammar[table];
	
	for (size_type i = 0; i != lattice.size(); ++ i)
	  if (! lattice[i].empty() && transducer.valid_span(i, i, 0)) {
	    edges.push_back(edge_type(semiring::traits<score_type>::one(), dot_type(table, transducer.root()), span_type(i, i)));
	    agenda_finishing[0].push_back(&edges.back());
	    //agenda_finishing.front().push_back(&edges.back());
	  }
      }

      // TODO: keep track of popped count for each span..
      // if we reached maximum, do not consider this span again...
      // how to control this...?
      
#if 0
      while (! agenda_finishing.front().empty()) {
	// explore traversals
	typename agenda_exploration_type::const_iterator aiter_end = agenda_exploration.end();
	for (typename agenda_exploration_type::const_iterator aiter = agenda_exploration.begin(); aiter != aiter_end; ++ aiter)
	  explore_traversal(*(*aiter), graph);
	agenda_exploration.clear();
	
	if (agenda_finishing.front().empty()) break;
	
	// we cannnot perform dynamic updates... this is TODO...
	agenda_finishing.front().make_heap();
	
	const edge_type* edge = agenda_finishing.front().top();
	agenda_finishing.front().pop();
	//std::cerr << "edge: span: " << edge->span << std::endl;
#if 0
	if (edge->is_passive())
	  std::cerr << "finished passive: " << edge->score
		    << " span: " << edge->span
		    << " lhs: " << edge->rule->rule->lhs << std::endl;
	else
	  std::cerr << "finished active: " << edge->score
		    << " span: " << edge->span
		    << " dot: " << edge->dot << std::endl;
#endif
	
	insert_hypergraph(*edge, lattice, graph);
	
#if 0
	if (edge->is_passive())
	  std::cerr << "graph size: " << graph.edges.size() << std::endl;
#endif
	
	if (edge->is_active()) {
	  if (pos_mode)
	    scan(*edge, lattice, __extract_terminal());
	  else
	    scan(*edge, lattice, __extract());
	  complete_active(*edge);
	} else {
	  complete_passive(*edge);
	  predict(*edge, lattice);
	}

	// check if we reached the goal!
	if (graph.is_valid()) break;
      }
#endif
      
#if 1
      const size_type num_popped_max = utils::bithack::max(lattice.size(), beam_size);
      
      for (;;) {
	// loop forever...
	
	for (size_type i = 0; i <= lattice.size(); ++ i) {

	  size_type num_popped = 0;
	  while (! agenda_finishing[i].empty() && num_popped != num_popped_max) {
	    // explore traversals
	    typename agenda_exploration_type::const_iterator aiter_end = agenda_exploration.end();
	    for (typename agenda_exploration_type::const_iterator aiter = agenda_exploration.begin(); aiter != aiter_end; ++ aiter)
	      explore_traversal(*(*aiter), graph);
	    agenda_exploration.clear();
	    
	    if (agenda_finishing[i].empty()) break;
	    
	    // we cannnot perform dynamic updates... this is TODO...
	    agenda_finishing[i].make_heap();
	    
	    const edge_type* edge = agenda_finishing[i].top();
	    agenda_finishing[i].pop();
	    num_popped += edge->is_passive();
	    
	    //std::cerr << "edge: span: " << edge->span << std::endl;
#if 0
	    if (edge->is_passive())
	      std::cerr << "finished passive: " << edge->score
			<< " span: " << edge->span
			<< " lhs: " << edge->rule->rule->lhs << std::endl;
	    else
	      std::cerr << "finished active: " << edge->score
			<< " span: " << edge->span
			<< " dot: " << edge->dot << std::endl;
#endif
	    
	    insert_hypergraph(*edge, lattice, graph);
	    	    
#if 0
	    if (edge->is_passive())
	      std::cerr << "graph size: " << graph.edges.size() << std::endl;
#endif
	    
	    if (edge->is_active()) {

	      if (pos_mode)
		scan(*edge, lattice, __extract_terminal());
	      else
		scan(*edge, lattice, __extract());
	      complete_active(*edge);
	    } else {
	      complete_passive(*edge);
	      predict(*edge, lattice);
	    }
	  }
	}
	
	// check if we reached the goal!
	if (graph.is_valid()) break;
      }
#endif
      
      if (graph.is_valid())
	graph.topologically_sort();
    }
    
  private:
    // TODO: scoring from passive edges should be accumulated...

    template <typename Extract>
    void scan(const edge_type& active, const lattice_type& lattice, Extract extractor)
    {
      if (active.span.last >= static_cast<int>(lattice.size())) return;

      //std::cerr << "scan active: " << active.span << ' ' << active.dot << std::endl;
      
      // scanning lattice...
      const transducer_type& transducer = grammar[active.dot.table];
      
      if (treebank && transducer.root() != active.dot.node) return;
      
      const lattice_type::arc_set_type& passive_arcs = lattice[active.span.last];
      
      const score_type score_prev = active.score / function(active.features);
      
      lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
      for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) {
	const symbol_type terminal = extractor(piter->label);

	//std::cerr << "terminal: " << terminal << std::endl;
	
	const span_type span_next(active.span.first, active.span.last + piter->distance);
	
	// handling of EPSILON rule...
	if (terminal == vocab_type::EPSILON) {
	  const dot_type& dot_next = active.dot;
	  const rule_candidate_ptr_set_type& rules = cands(active.dot.table, active.dot.node);

	  const feature_set_type features = active.features + piter->features;
	  const score_type       score_next = score_prev * function(features);
	  
	  // add passive edge... level is zero..
	  typename rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	  for (typename rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	    insert_edge(edge_type(score_next * function((*riter)->features), active, dot_next, span_next, *riter, features, active.attributes));
	  
	  // add active edge
	  insert_edge(edge_type(score_next, active, dot_next, span_next, features, active.attributes));
	  
	} else {
	  const transducer_type::id_type node = transducer.next(active.dot.node, terminal);
	  if (node == transducer.root()) continue;
	  
	  const dot_type dot_next(active.dot.table, node);
	  
	  //std::cerr << "dot next: " << dot_next << std::endl;

	  const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);
	  
	  const feature_set_type features = active.features + piter->features;
	  const score_type       score_next = score_prev * function(features);

	  //std::cerr << "rules size: " << rules.size() << std::endl;
	  
	  // add passive edge.. level is zero..
	  typename rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	  for (typename rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	    
	    //std::cerr << "scan: " << *((*riter)->rule)  << std::endl;
	    
	    insert_edge(edge_type(score_next * function((*riter)->features), active, dot_next, span_next, *riter, features, active.attributes));
	  }
	  
	  // add active edge
	  if (transducer.has_next(node))
	    insert_edge(edge_type(score_next, active, dot_next, span_next, features, active.attributes));
	}
      }
    }
    
    void predict(const edge_type& passive, const lattice_type& lattice)
    {
      if (passive.span.first == passive.span.last) return;
      //if (passive.span.level == 4) return; // max-unary chain

      //std::cerr << "predict passive: " << passive.span << std::endl;

      // from passive items, generate new actives...
      
      // extend root with passive items at [first, last)
      for (size_t table = 0; table != grammar.size(); ++ table) {
	const transducer_type& transducer = grammar[table];
	
	if (! transducer.valid_span(passive.span.first, passive.span.last, lattice.shortest_distance(passive.span.first, passive.span.last))) continue;
	if (! transducer.valid_span(passive.span.first, passive.span.first, 0)) continue;
	
	const transducer_type::id_type node = transducer.next(transducer.root(), passive.rule->rule->lhs);
	if (node == transducer.root()) continue;
	
	const dot_type dot_next(table, node);
	
	const rule_candidate_ptr_set_type& rules = cands(table, node);
	
	// add passive edge... this is definitedly unary-rule, thus level must be passive.level + 1!
	typename rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	for (typename rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	  const symbol_type& lhs = (*riter)->rule->lhs;
	  
	  const span_type span(passive.span.first, passive.span.last, utils::bithack::branch(lhs == goal, 0, passive.span.level + 1));
	  
	  //std::cerr << "predict: " << *((*riter)->rule)  << std::endl;
	  
	  insert_edge(edge_type(passive.score * function((*riter)->features), passive, dot_next, span, *riter));
	}
	
	// add active edge
	if (transducer.has_next(node))
	  insert_edge(edge_type(passive.score, passive, dot_next, passive.span));
      }
    }
    
    
    void complete_active(const edge_type& active)
    {
      if (active.span.first == active.span.last) return;

      //std::cerr << "complete active: " << active.span << std::endl;

      const transducer_type& transducer = grammar[active.dot.table];
      
      typename edge_ptr_set_type::const_iterator piter_end = edges_passive[active.span.last].end();
      for (typename edge_ptr_set_type::const_iterator piter = edges_passive[active.span.last].begin(); piter != piter_end; ++ piter) {
	const edge_type& passive = *(*piter);
	
	const transducer_type::id_type node = transducer.next(active.dot.node, passive.rule->rule->lhs);
	if (node == transducer.root()) continue;
#if 0
	std::cerr << "complete active: " << active.span << ' ' << active.dot
		  << " passive: " << passive.span << ' ' << passive.rule->rule->lhs
		  << std::endl;
#endif
	
	const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);

	const span_type  span_next(active.span.first, passive.span.last);
	const dot_type   dot_next(active.dot.table, node);
	const score_type score_next = active.score * passive.score;
	
	// add passive edge.. level is zero
	typename rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	for (typename rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	  //std::cerr << "complete active: " << *((*riter)->rule)  << std::endl;
	  
	  insert_edge(edge_type(score_next * function((*riter)->features), active, passive, dot_next, span_next, *riter, active.features, active.attributes));
	}
	
	// add active edge
	if (transducer.has_next(node))
	  insert_edge(edge_type(score_next, active, passive, dot_next, span_next, active.features, active.attributes));
      }
    }
    
    void complete_passive(const edge_type& passive)
    {
      if (passive.span.first == passive.span.last) return;
      
      //std::cerr << "complete passive: " << passive.span << std::endl;

      const symbol_type& lhs = passive.rule->rule->lhs;
      
      typename edge_ptr_set_type::const_iterator aiter_end = edges_active[passive.span.first].end();
      for (typename edge_ptr_set_type::const_iterator aiter = edges_active[passive.span.first].begin(); aiter != aiter_end; ++ aiter) {
	const edge_type& active = *(*aiter);
	
	const transducer_type& transducer = grammar[active.dot.table];
	
	const transducer_type::id_type node = transducer.next(active.dot.node, lhs);
	if (node == transducer.root()) continue;

#if 0
	std::cerr << "complete passive: " << passive.span << ' ' << passive.rule->rule->lhs
		  << " active: " << active.span << ' ' << active.dot
		  << std::endl;
#endif
	
	const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);
	
	const span_type  span_next(active.span.first, passive.span.last);
	const dot_type   dot_next(active.dot.table, node);
	const score_type score_next = active.score * passive.score;
	
	// add passive edge.. level is zero.
	typename rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	for (typename rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	  //std::cerr << "complete passive: " << *((*riter)->rule)  << std::endl;

	  insert_edge(edge_type(score_next * function((*riter)->features), active, passive, dot_next, span_next, *riter, active.features, active.attributes));
	}
	
	// add active edge
	if (transducer.has_next(node))
	  insert_edge(edge_type(score_next, active, passive, dot_next, span_next, active.features, active.attributes));
      }
    }
    
    void insert_edge(const edge_type& edge)
    {
      if (edge.passive && edge.active)
	if (! traversals.insert(traversal_type(edge.active, edge.passive, edge.is_active())).second)
	  return;
      
      edges.push_back(edge);
      
      agenda_exploration.push_back(&edges.back());
    }
    
    void explore_traversal(const edge_type& edge, hypergraph_type& graph)
    {
      if (edge.is_passive()) {
	typename discovered_passive_type::iterator diter = discovered_passive.find(&edge);
	if (diter == discovered_passive.end()) {
	  // newly discovered edge!
	  
	  // passive edge...
	  agenda_finishing[edge.span.size()].push_back(&edge);
	  //agenda_finishing.front().push_back(&edge);
	  
	  diter = discovered_passive.insert(std::make_pair(&edge, head_edge_set_type())).first;
	  diter->second.score = edge.score;
	  diter->second.edges.push_back(&edge);
	} else if (diter->second.head == hypergraph_type::invalid) {
	  diter->second.score = std::max(diter->second.score, edge.score);
	  diter->second.edges.push_back(&edge);
	} else {
	  // this will happen since we have already ignored potentially good edge!
	  diter->second.score = std::max(diter->second.score, edge.score);
	  insert_hypergraph(diter->second.head, edge, graph);
	}
      } else {
	if (discovered_active.insert(&edge).second) {
	  edges_active[edge.span.last].push_back(&edge);
	  
	  // active edge...
	  agenda_finishing[edge.span.size()].push_back(&edge);
	  //agenda_finishing.front().push_back(&edge);
	}
      }
    }
    
    void insert_hypergraph(const edge_type& edge, const lattice_type& lattice, hypergraph_type& graph)
    {
      // do not handle the special edge...
      if (edge.span.first == edge.span.last) return;

      if (! edge.is_passive()) return;
      
      edges_passive[edge.span.first].push_back(&edge);
      
      typename discovered_passive_type::iterator diter = discovered_passive.find(&edge);
      if (diter == discovered_passive.end()) {
	//std::cerr << "undiscovered passive?" << std::endl;
	return;
      }

      hypergraph_type::id_type head_id;
      if (diter->second.head == hypergraph_type::invalid) {
	if (edge.span.first == 0 && edge.span.last == static_cast<int>(lattice.size()) && edge.rule->rule->lhs == goal) {
	  if (! graph.is_valid())
	    graph.goal = graph.add_node().id;
	  head_id = graph.goal;
	} else
	  head_id = graph.add_node().id;
	
	diter->second.head = head_id;
      } else {
	//std::cerr << "already inserted?" << std::endl;
	head_id = diter->second.head;
      }
      
      typename edge_ptr_set_type::const_iterator eiter_end = diter->second.edges.end();
      for (typename edge_ptr_set_type::const_iterator eiter = diter->second.edges.begin(); eiter != eiter_end; ++ eiter)
	insert_hypergraph(head_id, *(*eiter), graph);
      
      const_cast<edge_ptr_set_type&>(diter->second.edges).clear();
    }

    std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails;

    void insert_hypergraph(const hypergraph_type::id_type& head_id, const edge_type& edge, hypergraph_type& graph)
    {
      tails.clear();
      
      const edge_type* curr = &edge;
      while (curr && grammar[curr->dot.table].root() != curr->dot.node) {
	if (curr->passive) {
	  typename discovered_passive_type::iterator diter = discovered_passive.find(curr->passive);
	  if (diter == discovered_passive.end())
	    throw std::runtime_error("no node?");
	    
	  tails.push_back(diter->second.head);
	}
	curr = curr->active;
      }
	
      std::reverse(tails.begin(), tails.end());
	
      hypergraph_type::edge_type& edge_new = graph.add_edge(tails.begin(), tails.end());
      edge_new.rule = edge.rule->rule;
      edge_new.features   = edge.features   + edge.rule->features;
      edge_new.attributes = edge.attributes + edge.rule->attributes;
      
      // assign metadata...
      edge_new.attributes[attr_span_first] = attribute_set_type::int_type(edge.span.first);
      edge_new.attributes[attr_span_last]  = attribute_set_type::int_type(edge.span.last);
      
      graph.connect_edge(edge_new.id, head_id);
    };
    
    const rule_candidate_ptr_set_type& cands(const size_type& table, const transducer_type::id_type& node)
    {
      typename rule_candidate_map_type::iterator riter = rule_tables[table].find(node);
      if (riter == rule_tables[table].end()) {
	riter = rule_tables[table].insert(std::make_pair(node, rule_candidate_ptr_set_type())).first;
	
	const transducer_type::rule_pair_set_type& rules = grammar[table].rules(node);
	
	riter->second.reserve(rules.size());
	
	transducer_type::rule_pair_set_type::const_iterator iter_begin = rules.begin();
	transducer_type::rule_pair_set_type::const_iterator iter_end   = rules.end();
	for (transducer_type::rule_pair_set_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter) {
	  rule_candidates.push_back(rule_candidate_type(function(iter->features),
							yield_source ? iter->source : iter->target,
							iter->features,
							iter->attributes));
	  
	  riter->second.push_back(&(rule_candidates.back()));
	}
      }
      return riter->second;
    }
    
  private:
    const symbol_type goal;
    const grammar_type& grammar;
    
    const function_type& function;
    
    const size_type  beam_size;
    const bool yield_source;
    const bool treebank;
    const bool pos_mode;
    const attribute_type attr_span_first;
    const attribute_type attr_span_last;
    
    edge_set_type edges;
    
    agenda_finishing_set_type agenda_finishing;
    agenda_exploration_type   agenda_exploration;
    
    traversal_set_type traversals;
    
    discovered_active_type  discovered_active;
    discovered_passive_type discovered_passive;
    
    edge_set_active_type    edges_active;
    edge_set_passive_type   edges_passive;
    
    rule_candidate_set_type   rule_candidates;
    rule_candidate_table_type rule_tables;
  };
  
  
  template <typename Function>
  inline
  void parse_agenda(const Symbol& goal, const Grammar& grammar, const Function& function, const Lattice& lattice, HyperGraph& graph, const int size, const bool yield_source=false, const bool treebank=false, const bool pos_mode=false)
  {
    ParseAgenda<typename Function::value_type, Function>(goal, grammar, function, size, yield_source, treebank, pos_mode)(lattice, graph);
  }

};

#endif
