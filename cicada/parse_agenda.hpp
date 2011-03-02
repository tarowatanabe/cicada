// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_AGENDA__HPP__
#define __CICADA__PARSE_AGENDA__HPP__ 1

#include <vector>
#include <algorithm>
#include <set>

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

#include <google/dense_hash_map>
#include <google/dense_hash_set>

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
    };
    
    typedef Dot dot_type;

    struct Span
    {
      int first;
      int last;
      
      Span() : first(0), last(0) {}
      Span(const int& __first, const int& __last) : first(__first), last(__last) {}
    };

    typedef Span span_type;

    struct Edge
    {
      typedef Edge edge_type;
      
      score_type  score;
      
      const edge_type* active;
      const edge_type* passive;
      
      dot_type                   dot;
      const rule_candidate_type* rule;
      feature_set_type           features;
      attribute_set_type         attributes;
      
      span_type   span;
      int         level;
      
    public:
      bool is_passive() const { return rule; }
      bool is_active() const { return ! rule; }
      
      bool is_scanned() const { return active && ! passive; }
      bool is_predicted() const { return ! active && passive; }
      bool is_completed() const { return active && passive; }
    };
    
    typedef Edge edge_type;
    
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

    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      if (lattice.empty()) return;
      
      // initialize by empty word... we do not support epsilon insertion as in the "Parsing and Hypergraphs"
      for (size_type table = 0; table != grammar.size(); ++ table) {
	const transducer_type& transducer = grammar[table];
	
	for (size_type i = 0; i != lattice.size(); ++ i)
	  if (! lattice[i].empty() && transducer.valid_span(i, i, 0))
	    insert_edge(edge_type(dot_type(table, transducer.root()), span_type(i, i)));
      }
	
      // main loop...
      while (! agenda_finishing.empty()) {
	// explore traversals
	agenda_type::const_iterator aiter_end = agenda_exploration.end();
	for (agenda_type::const_iterator aiter = agenda_exploration.begin(); aiter != aiter_end; ++ aiter)
	  explore_traversal(*(*aiter), source, target);
	agenda_exploration.clear();
	
	if (agenda_finishing.empty()) break;
	
	// we cannnot perform dynamic updates... this is TODO...
	agenda_finishing.make_heap();
	
	const edge_type* edge = agenda_finishing.top();
	agenda_finishing.pop();
	
	//
	// insert into graph... HOW?
	//
	
	if (edge->is_active()) {
	  scan(*edge);
	  complete_active(*edge);
	  predict(*edge);
	} else
	  complete_passive(*edge);
      }
    }
    
  private:
    void scan(const edge_tyep& active)
    {
      if (active.span.last >= lattice.size()) return;
      
      // scanning lattice...
      const transducer_type& transducer = grammar[active.dot.table];
      
      const lattice_type::arc_set_type& passive_arcs = lattice[active.span.last];
      
      lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
      for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) {
	const symbol_type& terminal = piter->label;
	
	const span_type span_next(active.span.first, active.span.last + piter->distance);
	
	// handling of EPSILON rule...
	if (terminal == vocab_type::EPSILON) {
	  const rule_candidate_ptr_set_type& rules = cands(active.dot.table, active.dot.node);
	  
	  // add passive edge... level is zero..
	  rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	  for (rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	    insert_edge(edge_type(active, dot_type(active.dot.table, active.dot.node), *riter, active.features + piter->features, active.attributes, span_next));
	  
	  // add active edge
	  insert_edge(edge_type(active, dot_type(active.dot.table, active.dot.node), active.features + piter->features, active.attributes, span_next));
	  
	} else {
	  const transducer_type::id_type node = transducer.next(active.dot.node, terminal);
	  if (node == transducer.root()) continue;
	  
	  const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);
	  
	  // add passive edge.. level is zero..
	  rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	  for (rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	    insert_edge(edge_type(active, dot_type(active.dot.table, node), *riter, active.features + piter->features, active.attributes, span_next));
	  
	  // add active edge
	  if (transducer.has_next(node))
	    insert_edge(edge_type(active, dot_type(active.dot.table, node), active.features + piter->featuers, active.attributes, span_next));
	}
      }
    }
    
    void predict(const edge_type& passive)
    {
      // from passive items, generate new actives...
      
      // extend root with passive items at [first, last)
      for (size_t table = 0; table != grammar.size(); ++ table) {
	const transducer_type& transducer = grammar[table];
	
	if (! transducer.valid_span(passive.span.first, passive.span.last, lattice.shortest_distance(passive.span.first, passive.span.last))) continue;
	if (! transducer.valid_span(passive.span.first, passive.span.first, 0)) continue;
	
	const transducer_type::id_type node = transducer.next(transducer.root(), passive.rule->rule->lhs);
	if (node == transducer.root()) continue;
	
	const span_type& span = passive.span;
	
	const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);
	
	// add passive edge... this is definitedly unary-rule, thus level must be passive.level + 1!
	rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	for (rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	  insert_edge(edge_type(passive, dot_type(table, node), *riter, span, passive.level + 1));
	
	// add active edge
	if (transducer.has_next(node))
	  insert_edge(edge_type(passive, dot_type(table, node), span));
      }
    }
    
    
    void complete_active(const edge_type& active)
    {
      const transducer_type& transducer = grammar[active.dot.table];

      const edge_type query(active.span.last, active.span.last);
      
      std::pair<edge_set_passive_type::const_iterator, edge_set_passive_type::const_iterator> result = edges_passive.equal_range(&query);
      for (edge_set_passive_type::const_iterator piter = result.first; piter != result.second; ++ piter) {
	const edge_type& passive = *(*piter);
	
	const transducer_type::id_type node = transducer.next(active.dot.node, non_terminal);
	if (node == transducer.root()) continue;
	
	const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);

	const span_type span_next(active.span.first, passive.span.last);
	
	// add passive edge.. level is zero
	rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	for (rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	  insert_edge(edge_type(active, passive, dot_type(active.dot.table, node), *riter, active.features, active.attributes, span_next));
	
	// add active edge
	if (transducer.has_next(node))
	  insert_edge(edge_type(active, passive, dot_type(active.dot.table, node), active.features, active.attributes, span_next));
      }
    }
    
    void complete_passive(const edge_type& passive)
    {
      const edge_type query(passive.span.first, passive.span.first);
      
      std::pair<edge_set_active_type::const_iterator, edge_set_active_type::const_iterator> result = edges_active.equal_range(&query);
      for (edge_set_active_type::const_iterator aiter = result.first; aiter != result.second; ++ aiter) {
	const edge_type& active = *(*aiter);
	
	const transducer_type& transducer = grammar[active.dot.table];
	
	const transducer_type::id_type node = transducer.next(active.dot.node, non_terminal);
	if (node == transducer.root()) continue;
	
	const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);
	
	const span_type span_next(active.span.first, passive.span.last);
	
	// add passive edge.. level is zero.
	rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	for (rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	  insert_edge(edge_type(active, passive, dot_type(active.dot.table, node), *riter, active.features, active.attributes, span_next));
	
	// add active edge
	if (transducer.has_next(node))
	  insert_edge(edge_type(active, passive, dot_type(active.dot.table, node), active.features, active.attributes, span_next));
      }
    }

    void insert_edge(const edge_type& edge)
    {
      
      
    }
    
    
    template <typename Tp>
    struct greater_ptr_score
    {
      bool operator()(const Tp* x, const Tp* y) const
      {
	return x->score > y->score;
      }
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
	
	std::sort(riter->second.begin(), riter->second.end(), greater_ptr_score<rule_candidate_type>());
      }
      return riter->second;
    }
    
  private:    
    
    rule_candidate_set_type   rule_candidates;
    rule_candidate_table_type rule_tables;
  };
};

#endif
