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
	
	insert_hypergraph(*edge, lattice, graph);
	
	// check if we reached the goal!
	if (graph.is_valid()) break;
	
	if (edge->is_active()) {
	  scan(*edge);
	  complete_active(*edge);
	  predict(*edge);
	} else
	  complete_passive(*edge);
      }
      
      if (graph.is_valid())
	graph.topologically_sort();
    }
    
  private:
    // TODO: scoring from passive edges should be accumulated...

    void scan(const edge_tyep& active)
    {
      if (active.span.last >= lattice.size()) return;
      
      // scanning lattice...
      const transducer_type& transducer = grammar[active.dot.table];
      
      const lattice_type::arc_set_type& passive_arcs = lattice[active.span.last];
      
      const score_type score_prev = active.score / function(active.features);
      
      lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
      for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) {
	const symbol_type& terminal = piter->label;
	
	const span_type span_next(active.span.first, active.span.last + piter->distance);
	
	// handling of EPSILON rule...
	if (terminal == vocab_type::EPSILON) {
	  const dot_type dot_next& dot_next = active.dot;
	  const rule_candidate_ptr_set_type& rules = cands(active.dot.table, active.dot.node);

	  const feature_set_type features = active.features + piter->features;
	  const score_type       score_next = score_prev * function(features);
	  
	  // add passive edge... level is zero..
	  rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	  for (rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	    insert_edge(edge_type(score_next * function((*riter)->features), active, dot_next, span_next, *riter, features, active.attributes));
	  
	  // add active edge
	  insert_edge(edge_type(score_next, active, dot_next, span_next, features, active.attributes));
	  
	} else {
	  const transducer_type::id_type node = transducer.next(active.dot.node, terminal);
	  if (node == transducer.root()) continue;

	  const dot_type dot_next(active.dot.table, node);
	  const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);
	  
	  const feature_set_type features = active.features + piter->features;
	  const score_type       score_next = score_prev * function(features);
	  
	  // add passive edge.. level is zero..
	  rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	  for (rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	    insert_edge(edge_type(score_next * function((*riter)->features), active, dot_next, span_next, *riter, features, active.attributes));
	  
	  // add active edge
	  if (transducer.has_next(node))
	    insert_edge(edge_type(score_next, active, dot_next, span_next, featuers, active.attributes));
	}
      }
    }
    
    void predict(const edge_type& passive)
    {
      if (passive.first == passive.last) return;

      // from passive items, generate new actives...
      
      // extend root with passive items at [first, last)
      for (size_t table = 0; table != grammar.size(); ++ table) {
	const transducer_type& transducer = grammar[table];
	
	if (! transducer.valid_span(passive.span.first, passive.span.last, lattice.shortest_distance(passive.span.first, passive.span.last))) continue;
	if (! transducer.valid_span(passive.span.first, passive.span.first, 0)) continue;
	
	const transducer_type::id_type node = transducer.next(transducer.root(), passive.rule->rule->lhs);
	if (node == transducer.root()) continue;
	
	const span_type& span = passive.span;
	const dot_type dot_next(table, node);
	
	const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);
	
	// add passive edge... this is definitedly unary-rule, thus level must be passive.level + 1!
	rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	for (rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	  const symbol_type& lhs = (*riter)->rule->lhs;
	  
	  insert_edge(edge_type(passive.score * function((*riter)->features), passive, dot_next, span, *riter, utils::bithack::branch(lhs == goal, 0, passive.level + 1)));
	}
	
	// add active edge
	if (transducer.has_next(node))
	  insert_edge(edge_type(semiring::traits<score_type>::one(), passive, dot_next, span));
      }
    }
    
    
    void complete_active(const edge_type& active)
    {
      if (passive.first == passive.last) return;

      const transducer_type& transducer = grammar[active.dot.table];

      const edge_type query(active.span.last, active.span.last);
      
      std::pair<edge_set_passive_type::const_iterator, edge_set_passive_type::const_iterator> result = edges_passive.equal_range(&query);
      for (edge_set_passive_type::const_iterator piter = result.first; piter != result.second; ++ piter) {
	const edge_type& passive = *(*piter);
	
	const transducer_type::id_type node = transducer.next(active.dot.node, non_terminal);
	if (node == transducer.root()) continue;
	
	const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);

	const span_type  span_next(active.span.first, passive.span.last);
	const dot_type   dot_next(active.dot.table, node);
	const score_type score_next = active.score * passive.score;
	
	// add passive edge.. level is zero
	rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	for (rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	  insert_edge(edge_type(score_next * function((*riter)->features), active, passive, dot_next, span_next, *riter, active.features, active.attributes));
	
	// add active edge
	if (transducer.has_next(node))
	  insert_edge(edge_type(score_next, active, passive, dot_next, span_next, active.features, active.attributes));
      }
    }
    
    void complete_passive(const edge_type& passive)
    {
      if (passive.first == passive.last) return;

      const edge_type query(passive.span.first, passive.span.first);
      
      std::pair<edge_set_active_type::const_iterator, edge_set_active_type::const_iterator> result = edges_active.equal_range(&query);
      for (edge_set_active_type::const_iterator aiter = result.first; aiter != result.second; ++ aiter) {
	const edge_type& active = *(*aiter);
	
	const transducer_type& transducer = grammar[active.dot.table];
	
	const transducer_type::id_type node = transducer.next(active.dot.node, non_terminal);
	if (node == transducer.root()) continue;
	
	const rule_candidate_ptr_set_type& rules = cands(active.dot.table, node);
	
	const span_type  span_next(active.span.first, passive.span.last);
	const dot_type   dot_next(active.dot.table, node);
	const score_type score_next = active.score * passive.score;
	
	// add passive edge.. level is zero.
	rule_candidate_ptr_set_type::const_iterator riter_end = rules.end();
	for (rule_candidate_ptr_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	  insert_edge(edge_type(score_next * function((*riter)->features)a, active, passive, dot_next, span_next, *riter, active.features, active.attributes));
	
	// add active edge
	if (transducer.has_next(node))
	  insert_edge(edge_type(score_next, active, passive, dot_next, span_next, active.features, active.attributes));
      }
    }
    
    void insert_edge(const edge_type& edge)
    {
      if (edge.passive && edge.active) {
	if (traversals.find(traversal_type(edge.active, edge.passive, edge.is_active())) != traversals.end())
	  return;
	else
	  traversals.insert(traversal_type(edge.active, edge.passive, edge.is_active()));
      }
      
      edges.push_back(edge);
      
      agenda_exploration.push_back(&edges.back());
    }
    
    void explore_traversal(const edge_type& edge)
    {
      if (edge.is_passive()) {
	discovered_type::const_iterator diter = discovered(edge.span.first, edges.span.last).find(symbol_level_type(edge.rule->rule->lhs, edge.level));
	
	if (diter == discovered.end()) {
	  // newly discovered edge!
	  agenda_finishing.push_back(&edges.back());
	  
	  discovered(edge.span.first, edges.span.last).insert(symbol_level_type(edge.rule->rule->lhs, edge.level));
	  
	  // add into hypergraph
	  graph_nodes[&edge].edges.push_back(&edge);
	} else {
	  // otherwise, updated scores in agenda_finishing...
	  edge_type& edge_discovered = const_cast<edge_type&>(*(*diter));
	  
	  edge_discovered.score = std::max(edge_discovered.score, edge.score);
	  
	  // add into hypergraph
	  graph_nodes[&edge_discovered].edges.push_back(&edge);
	}
      } else {
	if (edges_unique.find(&item) == edges_unique.end()) {
	  edges_unique.insert(&item);
	  
	  edges_active[edge.last].push_back(&edge);
	}
      }
    }
    
    void insert_hypergraph(const edge_type& edge, const lattice_type& lattice, hypergraph_type& graph)
    {
      // do not handle the special edge...
      if (edge.span.first == edge.span.last) return;
      
      edges_passive[edge.first].push_back(&edge);
      
      graph_node_set_type::iterator niter = graph_nodes.find(&edge);
      if (niter == graph_nodes.end()) continue;
      
      hypergraph_type::id_type head_id;
      if (edge.span.first == 0 && edge.span.last == lattie.size() && edge.rule->rule->lhs == goal) {
	if (! graph.is_valid())
	  graph.goal = graph.add_node().id;
	head_id = graph.goal;
      } else
	head_id = graph.add_node().id;
      
      niter->second.head = head_id;
      
      std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails;
      
      edge_ptr_set_type::const_iterator eiter_end = niter->second.edges.end();
      for (edge_ptr_set_type::const_iterator eiter = niter->second.edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge_antecedent = *(*eiter);
	
	tails.clear();
	
	const edge_type* curr = &edge_antecedent;
	while (curr && ! curr->is_predicted()) {
	  if (curr->passive) {
	    graph_node_set_type::iterator niter = graph_nodes.find(curr->passive);
	    if (niter == graph_nodes.end())
	      throw std::runtime_error("no node?");
	    
	    tails.push_back(niter->second.head);
	  }
	  curr = curr->active;
	}
	
	std::reverse(tails.begin(), tails.end());
	
	hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	edge_new.rule = edge_antecedent.rule->rule;
	edge_new.featuers   = edge_antecedent.features;
	edge_new.attributes = edge_antecedent.attributes;
	
	graph.connect_edge(edge_new.id, head_id);
      }
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
