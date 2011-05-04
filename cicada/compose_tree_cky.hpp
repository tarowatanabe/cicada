// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_TREE_CKY__HPP__
#define __CICADA__COMPOSE_TREE_CKY__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/sentence.hpp>
#include <cicada/vocab.hpp>
#include <cicada/tree_grammar.hpp>
#include <cicada/tree_transducer.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/lattice.hpp>
#include <cicada/hypergraph.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/bithack.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

//
// CFG parsing over lattice
//
// HOW TO HANDLE OOV????
//
// Insertion rule + glue rule...?
// [POS] ||| [x]
//  [x]  |||  x
// where x is a phrase...???
//

namespace cicada
{
  
  struct ComposeTreeCKY
  {
    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;
    
    typedef Sentence sentence_type;
    typedef Sentence phrase_type;

    typedef TreeGrammar    tree_grammar_type;
    typedef TreeTransducer tree_transducer_type;
    
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    
    typedef Lattice    lattice_type;
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef attribute_set_type::attribute_type attribute_type;
    
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    
    typedef tree_transducer_type::rule_type          tree_rule_type;
    typedef tree_transducer_type::rule_ptr_type      tree_rule_ptr_type;
    
    ComposeTreeCKY(const symbol_type& __goal, const tree_grammar_type& __tree_grammar, const grammar_type& __grammar, const bool __yield_source, const bool __unique_goal)
      : goal(__goal),
	tree_grammar(__tree_grammar), 
	grammar(__grammar),
	yield_source(__yield_source),
	unique_goal(__unique_goal),
	attr_span_first("span-first"),
	attr_span_last("span-last")
    {  
      goal_rule = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, goal.non_terminal())));
      
      closure.set_empty_key(symbol_type());
      closure_head.set_empty_key(symbol_type());
      closure_tail.set_empty_key(symbol_type());
    }
    
    struct ActiveTree
    {
      tree_transducer_type::id_type             node;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type   features;
      attribute_set_type attributes;
      
      ActiveTree(const tree_transducer_type::id_type& __node)
	: node(__node), tails(), features(), attributes() {}
      ActiveTree(const tree_transducer_type::id_type& __node,
		 const hypergraph_type::edge_type::node_set_type& __tails,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node), tails(__tails), features(__features), attributes(__attributes) {}
    };
    
    struct ActiveRule
    {
      transducer_type::id_type                  node;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type   features;
      attribute_set_type attributes;
      
      ActiveRule(const transducer_type::id_type& __node)
	: node(__node), tails(), features(), attributes() {}
      ActiveRule(const transducer_type::id_type& __node,
		 const hypergraph_type::edge_type::node_set_type& __tails,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node), tails(__tails), features(__features), attributes(__attributes) {}
    };

    typedef ActiveTree active_tree_type;
    typedef ActiveRule active_rule_type;
    
    typedef utils::chunk_vector<active_tree_type, 4096 / sizeof(active_tree_type), std::allocator<active_tree_type> > active_tree_set_type;
    typedef utils::chunk_vector<active_rule_type, 4096 / sizeof(active_rule_type), std::allocator<active_rule_type> > active_rule_set_type;
    
    typedef utils::chart<active_tree_set_type, std::allocator<active_tree_set_type> > active_tree_chart_type;
    typedef utils::chart<active_rule_set_type, std::allocator<active_rule_set_type> > active_rule_chart_type;
    
    typedef std::vector<active_tree_chart_type, std::allocator<active_tree_chart_type> > active_tree_chart_set_type;
    typedef std::vector<active_rule_chart_type, std::allocator<active_rule_chart_type> > active_rule_chart_set_type;
    
    typedef hypergraph_type::id_type passive_type;
    typedef std::vector<passive_type, std::allocator<passive_type> > passive_set_type;
    typedef utils::chart<passive_set_type, std::allocator<passive_set_type> > passive_chart_type;
    
    typedef std::pair<symbol_type, int> symbol_level_type;
    
    struct symbol_level_hash : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const symbol_level_type& x) const
      {
	return hasher_type::operator()(x.first, x.second);
      }
    };
    
    class NodeMap : public google::dense_hash_map<symbol_level_type, hypergraph_type::id_type, symbol_level_hash, std::equal_to<symbol_level_type> >
    {
    public:
      typedef google::dense_hash_map<symbol_level_type, hypergraph_type::id_type, symbol_level_hash, std::equal_to<symbol_level_type> > node_map_type;
      
    public:
      NodeMap() : node_map_type() { node_map_type::set_empty_key(symbol_level_type(symbol_type(), -1)); }
    };
    typedef NodeMap node_map_type;

    class NodeSet : public google::dense_hash_map<symbol_type, hypergraph_type::id_type, boost::hash<symbol_type>, std::equal_to<symbol_type> >
    {
    public:
      typedef google::dense_hash_map<symbol_type, hypergraph_type::id_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > node_set_type;
      
    public:
      NodeSet() : node_set_type() { node_set_type::set_empty_key(symbol_type()); }
    };
    typedef NodeSet node_set_type;
    
    typedef utils::chunk_vector<node_set_type, 4096 / sizeof(node_set_type), std::allocator<node_set_type> > node_graph_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
    
    typedef google::dense_hash_map<symbol_type, int, boost::hash<symbol_type>, std::equal_to<symbol_type> > closure_level_type;
    typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > closure_type;
    
    struct less_non_terminal
    {
      less_non_terminal(const non_terminal_set_type& __non_terminals) : non_terminals(__non_terminals) {}
      
      bool operator()(const hypergraph_type::id_type& x, const hypergraph_type::id_type& y) const
      {
	return non_terminals[x] < non_terminals[y] || (non_terminals[x] == non_terminals[y] && x < y);
      }
      
      const non_terminal_set_type& non_terminals;
    };
    
    void operator()(const lattice_type& lattice, hypergraph_type& graph)
    {
      graph.clear();
      
      if (lattice.empty())
	return;
      
      // initialize internal structure...
      actives_tree.clear();
      actives_rule.clear();
      passives.clear();
      
      node_map.clear();
      goal_id = hypergraph_type::invalid;
      non_terminals.clear();
      
      actives_tree.resize(tree_grammar.size(), active_tree_chart_type(lattice.size() + 1));
      actives_rule.resize(grammar.size(), active_rule_chart_type(lattice.size() + 1));
      
      passives.reserve(lattice.size() + 1);
      passives.resize(lattice.size() + 1);
      
      // initialize active chart
      for (size_t table = 0; table != tree_grammar.size(); ++ table) {
	const tree_transducer_type::id_type root = tree_grammar[table].root();
	
	for (size_t pos = 0; pos != lattice.size(); ++ pos)
	  actives_tree[table](pos, pos).push_back(active_tree_type(root));
      }
      
      for (size_t table = 0; table != grammar.size(); ++ table) {
	const transducer_type::id_type root = grammar[table].root();
	
	for (size_t pos = 0; pos != lattice.size(); ++ pos)
	  if (grammar[table].valid_span(pos, pos, 0))
	    actives_rule[table](pos, pos).push_back(active_rule_type(root));
      }
      
      for (size_t length = 1; length <= lattice.size(); ++ length)
	for (size_t first = 0; first + length <= lattice.size(); ++ first) {
	  const size_t last = first + length;
	  
	  node_map.clear();
	  
	  // tree-rules composed via CKY
	  for (size_t table = 0; table != tree_grammar.size(); ++ table) {
	    const tree_transducer_type& transducer = tree_grammar[table];
	    
	    // do we check by valid span...?
	    
	    // first, extend active items...
	    active_tree_set_type& cell = actives_tree[table](first, last);
	    for (size_t middle = first + 1; middle < last; ++ middle) {
	      const active_tree_set_type&  active_arcs  = actives_tree[table](first, middle);
	      const passive_set_type& passive_arcs = passives(middle, last);
	      
	      extend_actives(transducer, active_arcs, passive_arcs, cell);
	    }
	    
	    // then, advance by terminal(s) at lattice[last - 1];
	    {
	      const active_tree_set_type&       active_arcs  = actives_tree[table](first, last - 1);
	      const lattice_type::arc_set_type& passive_arcs = lattice[last - 1];
	      
	      active_tree_set_type::const_iterator aiter_begin = active_arcs.begin();
	      active_tree_set_type::const_iterator aiter_end   = active_arcs.end();
	      
	      if (aiter_begin != aiter_end) {
		lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
		for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) {
		  const symbol_type& terminal = piter->label;
		  
		  active_tree_set_type& cell = actives_tree[table](first, last - 1 + piter->distance);
		  
		  // handling of EPSILON rule...
		  if (terminal == vocab_type::EPSILON) {
		    for (active_tree_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
		      cell.push_back(active_type(aiter->node, aiter->tails, aiter->features + piter->features, aiter->attributes));
		  } else {
		    for (active_tree_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
		      const tree_transducer_type::id_type node = transducer.next(aiter->node, terminal);
		      if (node == transducer.root()) continue;
		      
		      cell.push_back(active_tree_type(node, aiter->tails, aiter->features + piter->features, aiter->attributes));
		    }
		  }
		}
	      }
	    }
	      
	    // complete active items if possible... The active items may be created from child span due to the
	    // lattice structure...
	    // apply rules on actives at [first, last)
	    
	    active_tree_set_type&  cell         = actives_tree[table](first, last);
	    passive_set_type&      passive_arcs = passives(first, last);
	    
	    active_tree_set_type::const_iterator citer_end = cell.end();
	    for (active_tree_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) {
	      const tree_transducer_type::rule_pair_set_type& rules = transducer.rules(citer->node);
	      
	      if (rules.empty()) continue;
	      
	      tree_transducer_type::rule_pair_set_type::const_iterator riter_begin = rules.begin();
	      tree_transducer_type::rule_pair_set_type::const_iterator riter_end   = rules.end();
	      
	      for (tree_transducer_type::rule_pair_set_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter) {
		const tree_rule_ptr_type& rule = (yield_source ? riter->source : riter->target);
		
		// apply tree-rule
		apply_rule(rule, riter->features + citer->features, riter->attributes + citer->attributes, citer->tails, 
			   passive_arcs, graph, first, last);
	      }
	    }
	  }
	  
	  // TODO: how to handle OOV???
	  // if we use grammar-insertion, rules cannot be instantiated...
	  // if we use generic POS symbol, then, we need to modify the symbol for the translational hypergraph...
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	      
	    // we will advance active spans, but constrained by transducer's valid span
	    // Here, we use only phrasal composition!
	    if (transducer.valid_span(first, last, lattice.shortest_distance(first, last))) {
	      
	      // first, extend active items...
	      active_rule_set_type& cell = actives_rule[table](first, last);
	      for (size_t middle = first + 1; middle < last; ++ middle) {
		const active_rule_set_type&  active_arcs  = actives_rule[table](first, middle);
		const passive_set_type& passive_arcs = passives(middle, last);
		
		extend_actives(transducer, active_arcs, passive_arcs, cell);
	      }
	      
	      // then, advance by terminal(s) at lattice[last - 1];
	      const active_set_type&  active_arcs  = actives_rule[table](first, last - 1);
	      const lattice_type::arc_set_type& passive_arcs = lattice[last - 1];
	      
	      active_set_type::const_iterator aiter_begin = active_arcs.begin();
	      active_set_type::const_iterator aiter_end = active_arcs.end();
	      
	      if (aiter_begin != aiter_end) {
		lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
		for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) {
		  const symbol_type& terminal = piter->label;
		  
		  active_rule_set_type& cell = actives_rule[table](first, last - 1 + piter->distance);
		  
		  // handling of EPSILON rule...
		  if (terminal == vocab_type::EPSILON) {
		    for (active_rule_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
		      cell.push_back(active_rule_type(aiter->node, aiter->features + piter->features, aiter->attributes));
		  } else {
		    for (active_rule_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
		      const transducer_type::id_type node = transducer.next(aiter->node, terminal);
		      if (node == transducer.root()) continue;
		      
		      cell.push_back(active_rule_type(node, aiter->features + piter->features, aiter->attributes));
		    }
		  }
		}
	      }
	    }
	    
	    // complete...active items...
	    // when we found rules, we try all possible combination of lhs already used for this span!
	    
	    active_rule_set_type&  cell         = actives_rule[table](first, last);
	    passive_set_type&      passive_arcs = passives(first, last);

	    // when passive_arcs is empty, we have no span...! How to handle this case...???
	    
	    active_rule_set_type::const_iterator citer_end = cell.end();
	    for (active_rule_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) {
	      const transducer_type::rule_pair_set_type& rules = transducer.rules(citer->node);
	      
	      if (rules.empty()) continue;
	      
	      transducer_type::rule_pair_set_type::const_iterator riter_begin = rules.begin();
	      transducer_type::rule_pair_set_type::const_iterator riter_end   = rules.end();
	      
	      for (transducer_type::rule_pair_set_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter) {
		const rule_ptr_type& rule = (yield_source ? riter->source : riter->target);
		
		// here, we will try match with all the non-terminals for [first, last), ignoring this rule's lhs
		passive_set_type::const_iterator piter_begin = passives.begin();
		passive_set_type::const_iterator piter_end   = passives.end();
		
		for (passive_set_type::const_iterator piter = piter_begin; piter != piter_end; ++ piter) {
		  const symbol_type& lhs = non_terminals[*piter];
		  
		  edge_type& edge = graph.add_edge();
		  edge.rule = rule_type::create(rule_type(lhs, rule->rhs));
		  edge.features   = riter->features;
		  edge.attributes = riter->attributes;
		  
		  graph.connect_edge(edge.id, *piter);
		}
	      }
	    }
	  }
	  
	  // handle unary rules...
	  // TODO: handle unary rules both for tree-grammar and grammar!!!!
	  if (! passives(first, last).empty()) {
	    //std::cerr << "closure from passives: " << passives(first, last).size() << std::endl;
	    
	    passive_set_type& passive_arcs = passives(first, last);
	    
	    size_t passive_first = 0;
	    
	    // initialize closure..
	    closure.clear();
	    passive_set_type::const_iterator piter_end = passive_arcs.end();
	    for (passive_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter)
	      closure[non_terminals[*piter]] = 0;

	    edge_type::node_set_type tails(1);
	    
	    int unary_loop = 0;
	    for (;;) {
	      const size_t passive_size = passive_arcs.size();
	      const size_t closure_size = closure.size();
	      
	      closure_head.clear();
	      closure_tail.clear();
	      
	      for (size_t table = 0; table != tree_grammar.size(); ++ table) {
		const tree_transducer_type& transducer = tree_grammar[table];
		
		for (size_t p = passive_first; p != passive_size; ++ p) {
		  const symbol_type non_terminal = non_terminals[passive_arcs[p]];
		  
		  const tree_transducer_type::id_type node = transducer.next(transducer.root(), non_terminal);
		  if (node == transducer.root()) continue;
		  
		  const tree_transducer_type::rule_pair_set_type& rules = transducer.rules(node);
		  
		  if (rules.empty()) continue;
		  
		  closure_tail.insert(non_terminal);
		  
		  tree_transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
		  for (tree_transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
		    const tree_rule_ptr_type& rule = (yield_source ? riter->source : riter->target);
		    const symbol_type& lhs = rule->label;
		    
		    closure_level_type::const_iterator citer = closure.find(lhs);
		    const int level = (citer != closure.end() ? citer->second : 0);
		    
		    closure_head.insert(lhs);
		    
		    tails.front() = passive_arcs[p];
		    
		    // apply rule...
		    apply_rule(rule, riter->features, riter->attributes, tails, 
			       passive_arcs, graph, first, last, level + 1);
		    
		  }
		}
	      }

	      for (size_t table = 0; table != grammar.size(); ++ table) {
		const transducer_type& transducer = grammar[table];
		
		if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
		
		for (size_t p = passive_first; p != passive_size; ++ p) {
		  const symbol_type non_terminal = non_terminals[passive_arcs[p]];
		  
		  const transducer_type::id_type node = transducer.next(transducer.root(), non_terminal);
		  if (node == transducer.root()) continue;
		  
		  const transducer_type::rule_pair_set_type& rules = transducer.rules(node);
		  
		  if (rules.empty()) continue;
		  
		  // passive_arcs "MAY" be modified!
		  
		  closure_tail.insert(non_terminal);
		  
		  transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
		  for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
		    const rule_ptr_type& rule = (yield_source ? riter->source : riter->target);
		    const symbol_type& lhs = rule->lhs;
		    
		    closure_level_type::const_iterator citer = closure.find(lhs);
		    const int level = (citer != closure.end() ? citer->second : 0);
		    
		    closure_head.insert(lhs);
		    
		    tails.front() = passive_arcs[p];
		    
		    apply_rule(rule, riter->features, riter->attributes, tails,
			       passive_arcs, graph, first, last, level + 1);
		  }
		}
	      }
	      
	      if (passive_size == passive_arcs.size()) break;
	      
	      passive_first = passive_size;
	      
	      // we use level-one, that is the label assigned for new-lhs!
	      closure_type::const_iterator hiter_end = closure_head.end();
	      for (closure_type::const_iterator hiter = closure_head.begin(); hiter != hiter_end; ++ hiter)
		closure.insert(std::make_pair(*hiter, 1));
	      
	      // increment non-terminal level when used as tails...
	      closure_type::const_iterator titer_end = closure_tail.end();
	      for (closure_type::const_iterator titer = closure_tail.begin(); titer != titer_end; ++ titer)
		++ closure[*titer];
	      
	      if (closure_size != closure.size())
		unary_loop = 0;
	      else
		++ unary_loop;
	      
	      // 4 iterations
	      if (unary_loop == 4) break;
	    }
	  }
	  
	  // sort passives at passives(first, last) wrt non-terminal label in non_terminals
	  std::sort(passives(first, last).begin(), passives(first, last).end(), less_non_terminal(non_terminals));
	  
	  // extend root with passive items at [first, last)
	  // we need to do this for simple transducers, also...
	  for (size_t table = 0; table != tree_grammar.size(); ++ table) {
	    const tree_transducer_type& transducer = tree_grammar[table];
	    
	    const active_tree_set_type& active_arcs  = actives_tree[table](first, first);
	    const passive_set_type&     passive_arcs = passives(first, last);
	    
	    active_tree_set_type& cell = actives_tree[table](first, last);
	    
	    extend_actives(transducer, active_arcs, passive_arcs, cell);
	  }
	  
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	    
	    const active_rule_set_type& active_arcs  = actives_rule[table](first, first);
	    const passive_set_type&     passive_arcs = passives(first, last);
	    
	    active_rule_set_type& cell = actives_rule[table](first, last);
	    
	    extend_actives(transducer, active_arcs, passive_arcs, cell);
	  }
	  
	  // final patch work...
	  //
	  // From the pool of actives in actives_rule[table](first, last)
	  // find phrasal rule, and try match with existing forest if not already used!
	  // we will be flexible in adjusting the lhs symbol so that we can match it with non-temrinals in the forest.
	  
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	    
	    active_rule_set_type&  cell         = actives_rule[table](first, last);
	    passive_set_type&      passive_arcs = passives(first, last);
	    
	    // we will collect only phrasal rules...
	    active_rule_set_type::const_iterator citer_end = cell.end();
	    for (active_rule_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) 
	      if (citer->tails.empty()) {
		const transducer_type::rule_pair_set_type& rules = transducer.rules(citer->node);
		
		if (rules.empty()) continue;
		
		transducer_type::rule_pair_set_type::const_iterator riter_begin = rules.begin();
		transducer_type::rule_pair_set_type::const_iterator riter_end   = rules.end();
		
		for (transducer_type::rule_pair_set_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter) {
		  const rule_ptr_type& rule = (yield_source ? riter->source : riter->target);
		  
		  if (node_map_local.find(rule->lhs) != node_map_local.end()) continue;
		  
		  node_map_type::const_iterator liter_end = node_map_local.end();
		  for (node_map_type::const_iterator liter = node_map_local.begin(); liter != liter_end; ++ liter) 
		    if (liter->first != rule->lhs) {
		      const symbol_type& lhs = liter->first;
		      
		      // create new rule
		      edge_type& edge_new = graph.add_edge();
		      edge_new.rule = rule_type::create(rule_type(lhs, rule->rhs));
		      edge_new.features   = riter->features;
		      edge_new.attributes = riter->attributes;
		      
		      edge.attributes[attr_span_first] = attribute_set_type::int_type(first);
		      edge.attributes[attr_span_last]  = attribute_set_type::int_type(last);
		    }
		}
	      }
	  }
	}
      
      // final...
      
      // we will clear node map so that we will always create new node..
      
      node_map.clear();
      
      if (unique_goal) {
	passive_set_type& passive_arcs = passives(0, lattice.size());
	for (size_t p = 0; p != passive_arcs.size(); ++ p)
	  if (non_terminals[passive_arcs[p]] == goal) {
	    if (graph.is_valid())
	      throw std::runtime_error("multiple goal? " + boost::lexical_cast<std::string>(graph.goal) + " " + boost::lexical_cast<std::string>(passive_arcs[p]));
	    
	    graph.goal = passive_arcs[p];
	  }
	
      } else {
	edge_type::node_set_type tails(1);
	passive_set_type& passive_arcs = passives(0, lattice.size());
	for (size_t p = 0; p != passive_arcs.size(); ++ p)
	  if (non_terminals[passive_arcs[p]] == goal) {
	    tails.front() = passive_arcs[p];
	    
	    apply_rule(goal_rule, feature_set_type(), attribute_set_type(), tails, 
		       passive_arcs, graph, 0, lattice.size(), 0, true);
	  }
      }
    }
    
    template <typename Transducer, typename Actives>
    bool extend_actives(const Transducer& transducer,
			const Actives& actives, 
			const passive_set_type& passives,
			Actives& cell)
    {
      typename Actives::const_iterator aiter_begin = actives.begin();
      typename Actives::const_iterator aiter_end   = actives.end();
      
      passive_set_type::const_iterator piter_begin = passives.begin();
      passive_set_type::const_iterator piter_end   = passives.end();
      
      bool found = false;
      
      if (piter_begin != piter_end)
	for (typename Actives::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  if (transducer.has_next(aiter->node)) {
	    symbol_type label;
	    typename Transducer::id_type node = transducer.root();
	    
	    hypergraph_type::edge_type::node_set_type tails(aiter->tails.size() + 1);
	    std::copy(aiter->tails.begin(), aiter->tails.end(), tails.begin());
	    
	    for (passive_set_type::const_iterator piter = piter_begin; piter != piter_end; ++ piter) {
	      const symbol_type& non_terminal = non_terminals[*piter];
	      
	      if (label != non_terminal) {
		node = transducer.next(aiter->node, non_terminal);
		label = non_terminal;
	      }
	      if (node == transducer.root()) continue;
	      
	      tails.back() = *piter;
	      cell.push_back(typename Actives::value_type(node, tails, aiter->features, aiter->attributes));
	      
	      found = true;
	    }
	  }
      
      return found;
    }
    
    
    void apply_rule(const symbol_type& lhs,
		    const rule_ptr_type& rule,
		    const feature_set_type& features,
		    const attribute_set_type& attributes,
		    const hypergraph_type::edge_type::node_set_type& frontier,
		    passive_set_type& passives,
		    hypergraph_type& graph,
		    const int lattice_first,
		    const int lattice_last,
		    const int level = 0,
		    const bool is_goal = false)
    {
      hypergraph_type::edge_type& edge = graph.add_edge(frontier);
      edge.rule = rule;
      edge.features   = features;
      edge.attributes = attributes;
      
      // assign metadata...
      edge.attributes[attr_span_first] = attribute_set_type::int_type(lattice_first);
      edge.attributes[attr_span_last]  = attribute_set_type::int_type(lattice_last);

      if (is_goal) {
	if (! graph.is_valid()) {
	  graph.goal = graph.add_node().id;
	  non_terminals.push_back(rule->lhs);
	}
	
	graph.connect_edge(edge.id, graph.goal);
      } else {
	const int cat_level = utils::bithack::branch(unique_goal && rule->lhs == goal, 0, level);
	
	std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(rule->lhs, cat_level), 0));
	if (result.second) {
	  hypergraph_type::node_type& node = graph.add_node();
	  non_terminals.push_back(rule->lhs);
	  passives.push_back(node.id);
	  result.first->second = node.id;
	}
	
	graph.connect_edge(edge.id, result.first->second);
      }
    }
    
    void apply_rule(const symbol_type& lhs,
		    const tree_rule_type& rule,
		    const feature_set_type& features,
		    const attribute_set_type& attributes,
		    const hypergraph_type::edge_type::node_set_type& frontier,
		    passive_set_type& passives,
		    hypergraph_type& graph,
		    const int lattice_first,
		    const int lattice_last,
		    const int level = 0,
		    const bool is_goal = false)
    {
      hypergrap_type::id_type root_id = hypergraph_type::invalid;
      if (is_goal) {
	if (! graph.is_valid()) {
	  graph.goal = graph.add_node().id;
	  non_termnals.push_back(rule->lhs);
	}
	root_id = graph.goal;
      } else {
	const int cat_level = utils::bithack::branch(unique_goal && rule->lhs == goal, 0, level);
	
	std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(rule->lhs, cat_level), 0));
	if (result.second) {
	  hypergraph_type::node_type& node = graph.add_node();
	  non_terminals.push_back(rule->lhs);
	  passives.push_back(node.id);
	  result.first->second = node.id;
	}
	root_id = result.first->second;
      }
      
      int non_terminal_pos = 0;
      
      const hypergraph_type::id_type edge_id = construct_graph(rule, root_id, frontier, graph, non_terminal_pos);
      
      graph.edges[edge_id].features   = features;
      graph.edges[edge_id].attributes = attributes;
      
      // assign metadata only for the root edge...?????
      graph.edges[edge_id].attributes[attr_span_first] = attribute_set_type::int_type(lattice_first);
      graph.edges[edge_id].attributes[attr_span_last]  = attribute_set_type::int_type(lattice_last);
    }
    
    hypergraph_type::id_type construct_graph(const tree_rule_type& rule,
					     const hypergraph_type::id_type root,
					     const hypergraph_type::edge_type::node_set_type& frontiers,
					     hypergraph_type& graph,
					     int& non_terminal_pos)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > rhs_type;
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails_type;
      
      rhs_type rhs;
      tails_type tails;
      
      int pos = 1;
      tree_rule_type::const_iterator aiter_end = rule.end();
      for (tree_rule_type::const_iterator aiter = rule.begin(); aiter != aiter_end; ++ aiter)
	if (aiter->label.is_non_terminal()) {
	  if (aiter->antecedents.empty()) {
	    const int __non_terminal_index = aiter->label.non_terminal_index();
	    const int non_terminal_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	    ++ non_terminal_pos;
	    
	    if (non_terminal_index >= static_cast<int>(frontiers.size()))
	      throw std::runtime_error("non-terminal index exceeds frontier size");
	    
	    tails.push_back(frontiers[non_terminal_index]);
	  } else
	    tails.push_back(graph.add_node().id); // tree-rule private node!
	  
	  rhs.push_back(aiter->label.non_terminal(pos));
	  
	  ++ pos;
	} else
	  rhs.push_back(aiter->label);
      
      const hypergraph_type::id_type edge_id = graph_out.add_edge(tails.begin(), tails.end()).id;
      graph_out.edges[edge_id].rule = rule_type::create(rule_type(rule.label, rhs.begin(), rhs.end()));
      
      graph_out.connect_edge(edge_id, root);
      
      // if we have antecedents traverse and construct
      tails_type::const_iterator titer = tails.begin();
      for (tree_rule_type::const_iterator aiter = rule.begin(); aiter != aiter_end; ++ aiter)
	if (aiter->label.is_non_terminal() && ! aiter->antecedents.empty()) {
	  construct_graph(*aiter, *titer, frontiers, graph, non_terminal_pos);
	  ++ titer;
	}
      
      return edge_id;
    }
    
  private:
    const symbol_type goal;
    const tree_grammar_type& tree_grammar;
    const grammar_type& grammar;
    const bool yield_source;
    const bool unique_goal;
    
    attribute_type attr_span_first;
    attribute_type attr_span_last;
    
    rule_ptr_type goal_rule;
    
    active_tree_chart_set_type actives_tree;
    active_rule_chart_set_type actives_rule;
    passive_chart_type         passives;
    
    closure_level_type    closure;
    closure_type          closure_head;
    closure_type          closure_tail;
    
    node_map_type           node_map;
    node_graph_type         graph_source;
    hypegraph_type::id_type graph_goal;
    non_terminal_set_type   non_terminals;
  };
  
  
  inline
  void compose_tree_cky(const Symbol& goal, const TreeGrammar& tree_grammar, const Grammar& grammar, const Lattice& lattice, HyperGraph& graph, const bool yield_source=false, const bool unique_goal=false)
  {
    ComposeTreeCKY __composer(goal, tree_grammar, grammar, yield_source, unique_goal);
    __composer(lattice, graph);
  }
};

#endif
