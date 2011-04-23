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
    
    typedef tree_transducer_type::rule_pair_set_type tree_rule_pair_set_type;
    typedef tree_transducer_type::rule_pair_type     tree_rule_pair_type;
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
      
      node_map.set_empty_key(symbol_level_type());
      closure.set_empty_key(symbol_type());
      closure_head.set_empty_key(symbol_type());
      closure_tail.set_empty_key(symbol_type());
    }
    
    struct ActiveTree
    {
      
      tree_transducer_type::id_type             node;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type                          features;
      attribute_set_type                        attributes;
    };
    
    struct ActivePhrase
    {
      transducer_type::id_type node;
      feature_set_type         features;
      attribute_set_type       attributes;
    };

    typedef ActiveTree   active_tree_type;
    typedef ActivePhrase active_phrase_type;
    
    typedef utils::chunk_vector<active_tree_type, 4096 / sizeof(active_tree_type), std::allocator<active_tree_type> >       active_tree_set_type;
    typedef utils::chunk_vector<active_phrase_type, 4096 / sizeof(active_phrase_type), std::allocator<active_phrase_type> > active_phrase_set_type;
    
    typedef utils::chart<active_tree_set_type, std::allocator<active_tree_set_type> >     active_tree_chart_type;
    typedef utils::chart<active_phrase_set_type, std::allocator<active_phrase_set_type> > active_phrase_chart_type;
    
    typedef std::vector<active_tree_chart_type, std::allocator<active_tree_chart_type> >     active_tree_chart_set_type;
    typedef std::vector<active_phrase_chart_type, std::allocator<active_phrase_chart_type> > active_phrase_chart_set_type;
    
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
    
    typedef google::dense_hash_map<symbol_level_type, hypergraph_type::id_type, symbol_level_hash, std::equal_to<symbol_level_type> > node_map_type;
    
    typedef google::dense_hash_map<symbol_type, int, boost::hash<symbol_type>, std::equal_to<symbol_type> > closure_level_type;
    typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > closure_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
    
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
      actives_phrase.clear();
      passives.clear();
      non_terminals.clear();
      
      actives_tree.resize(tree_grammar.size(), active_tree_chart_type(lattice.size() + 1));
      actives_phrase.resize(grammar.size(), active_phrase_chart_type(lattice.size() + 1));
      passives.resize(lattice.size() + 1);
      
      // initialize active chart
      for (size_t table = 0; table != tree_grammar.size(); ++ table) {
	const transducer_type::id_type root = tree_grammar[table].root();
	
	for (size_t pos = 0; pos != lattice.size(); ++ pos)
	  actives_tree[table](pos, pos).push_back(active_tree_type(root));
      }
      
      for (size_t table = 0; table != grammar.size(); ++ table) {
	const transducer_type::id_type root = grammar[table].root();
	
	for (size_t pos = 0; pos != lattice.size(); ++ pos)
	  if (grammar[table].valid_span(pos, pos, 0))
	    actives_phrase[table](pos, pos).push_back(active_phrase_type(root));
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
		// apply tree-rule
		
	      }
	    }
	  }
	  
	  // handle unary rules...
	  if (! passives(first, last).empty()) {
	    
	    
	  }
	  
	  
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	    
	    // we will advance active spans, but constrained by transducer's valid span
	    // Here, we use only phrasal composition!
	    if (transducer.valid_span(first, last, lattice.shortest_distance(first, last))) {
	      
	      // advance by terminal(s) at lattice[last - 1];
	      const active_set_type&  active_arcs  = actives_phrase[table](first, last - 1);
	      const lattice_type::arc_set_type& passive_arcs = lattice[last - 1];
	      
	      active_set_type::const_iterator aiter_begin = active_arcs.begin();
	      active_set_type::const_iterator aiter_end = active_arcs.end();
	      
	      if (aiter_begin != aiter_end) {
		lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
		for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) {
		  const symbol_type& terminal = piter->label;
		  
		  active_phrase_set_type& cell = actives_phrase[table](first, last - 1 + piter->distance);
		  
		  // handling of EPSILON rule...
		  if (terminal == vocab_type::EPSILON) {
		    for (active_phrase_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
		      cell.push_back(active_phrase_type(aiter->node, aiter->features + piter->features, aiter->attributes));
		  } else {
		    for (active_phrase_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
		      const transducer_type::id_type node = transducer.next(aiter->node, terminal);
		      if (node == transducer.root()) continue;
		      
		      cell.push_back(active_phrase_type(node, aiter->features + piter->features, aiter->attributes));
		    }
		  }
		}
	      }
	    }
	    
	    // complete...active items...
	    // when we found phrases, we try all possible combination of lhs already used for this span!
	    
	    active_phrase_set_type&  cell         = actives_phrase[table](first, last);
	    passive_set_type&        passive_arcs = passives(first, last);
	    
	    active_phrase_set_type::const_iterator citer_end = cell.end();
	    for (active_phrase_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) {
	      const transducer_type::rule_pair_set_type& rules = transducer.rules(citer->node);
	      
	      if (rules.empty()) continue;
	      
	      transducer_type::rule_pair_set_type::const_iterator riter_begin = rules.begin();
	      transducer_type::rule_pair_set_type::const_iterator riter_end   = rules.end();
	      
	      for (transducer_type::rule_pair_set_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter) {
		const rule_ptr_type& rule = (yield_source ? riter->source : riter->target);
		
		// here, we will try match with all the non-terminals for [first, last), ignoring this rule's lhs
		
		
	      }
	    }
	  }
	  
	  // sort passives at passives(first, last) wrt non-terminal label in non_terminals
	  std::sort(passives(first, last).begin(), passives(first, last).end(), less_non_terminal(non_terminals));
	  
	  // extend root with passive items at [first, last)
	  for (size_t table = 0; table != tree_grammar.size(); ++ table) {
	    const tree_transducer_type& transducer = tree_grammar[table];
	    
	    const active_tree_set_type& active_arcs  = actives_tree[table](first, first);
	    const passive_set_type&     passive_arcs = passives(first, last);
	    
	    active_tree_set_type& cell = actives_tree[table](first, last);
	    
	    extend_actives(transducer, active_arcs, passive_arcs, cell);
	  }
	}
      
      // final...
    }
    
    const symbol_type goal;
    const tree_grammar_type& tree_grammar;
    const grammar_type& grammar;
    const bool yield_source;
    const bool unique_goal;
    
    attribute_type attr_span_first;
    attribute_type attr_span_last;
    
    rule_ptr_type goal_rule;
    
    active_tree_chart_set_type   actives_tree;
    active_phrase_chart_set_type actives_phrase;
    passive_chart_type           passives;

    node_map_type         node_map;
    closure_level_type    closure;
    closure_type          closure_head;
    closure_type          closure_tail;
    non_terminal_set_type non_terminals;
  };
  
  
  inline
  void compose_tree_cky(const Symbol& goal, const TreeGrammar& tree_grammar, const Grammar& grammar, const Lattice& lattice, HyperGraph& graph, const bool yield_source=false, const bool unique_goal=false)
  {
    ComposeTreeCKY __composer(goal, tree_grammar, grammar, yield_source, unique_goal);
    __composer(lattice, graph);
  }
};

#endif
