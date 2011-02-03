// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_CKY__HPP__
#define __CICADA__COMPOSE_CKY__HPP__ 1

#include <deque>
#include <vector>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/hashmurmur.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

namespace cicada
{
  
  struct ComposeCKY
  {
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

    
    ComposeCKY(const symbol_type& __goal, const grammar_type& __grammar, const bool __yield_source=false)
      : goal(__goal), grammar(__grammar), yield_source(__yield_source),
	attr_span_first("span-first"),
	attr_span_last("span-last")
    {
      goal_rule = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, goal.non_terminal(1))));
      
      node_map.set_empty_key(symbol_index_type());
      closure.set_empty_key(symbol_type());
    }
    
    struct ActiveItem
    {
      ActiveItem(const transducer_type::id_type& __node,
		 const hypergraph_type::edge_type::node_set_type __tails,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node),
	  tails(__tails),
	  features(__features),
	  attributes(__attributes) {}
      ActiveItem(const transducer_type::id_type& __node,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node),
	  tails(),
	  features(__features),
	  attributes(__attributes) {}
      ActiveItem(const transducer_type::id_type& __node,
		 const hypergraph_type::edge_type::node_set_type __tails,
		 const feature_set_type& __features)
	: node(__node),
	  tails(__tails),
	  features(__features),
	  attributes() {}
      ActiveItem(const transducer_type::id_type& __node,
		 const feature_set_type& __features)
	: node(__node),
	  tails(),
	  features(__features),
	  attributes() {}
      ActiveItem(const transducer_type::id_type& __node)
	: node(__node),
	  tails(),
	  features(),
	  attributes() {}
      
      transducer_type::id_type                  node;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type                          features;
      attribute_set_type                        attributes;
    };
    
    typedef ActiveItem active_type;
    typedef utils::chunk_vector<active_type, 4096 / sizeof(active_type), std::allocator<active_type> > active_set_type;
    typedef utils::chart<active_set_type, std::allocator<active_set_type> > active_chart_type;
    typedef std::vector<active_chart_type, std::allocator<active_chart_type> > active_chart_set_type;

    typedef hypergraph_type::id_type passive_type;
    typedef std::vector<passive_type, std::allocator<passive_type> > passive_set_type;
    typedef utils::chart<passive_set_type, std::allocator<passive_set_type> > passive_chart_type;

    typedef std::pair<symbol_type, int> symbol_index_type;
    typedef google::dense_hash_map<symbol_index_type, hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<symbol_index_type> > node_map_type;
    
    typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > closure_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
    
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      if (lattice.empty())
	return;
      
      // initialize internal structure...
      actives.clear();
      passives.clear();
      non_terminals.clear();
      
      actives.resize(grammar.size(), active_chart_type(lattice.size() + 1));
      passives.resize(lattice.size() + 1);
      
      // initialize active chart
      for (size_t table = 0; table != grammar.size(); ++ table) {
	const transducer_type::id_type root = grammar[table].root();
	
	for (size_t pos = 0; pos != lattice.size(); ++ pos)
	  if (grammar[table].valid_span(pos, pos, 0))
	    actives[table](pos, pos).push_back(active_type(root));
      }

      
      for (size_t length = 1; length <= lattice.size(); ++ length)
	for (size_t first = 0; first + length <= lattice.size(); ++ first) {
	  const size_t last = first + length;
	  
	  node_map.clear();
	  
	  //std::cerr << "span: " << first << ".." << last << " distance: " << lattice.shortest_distance(first, last) << std::endl;
	  
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	    
	    // we will advance active spans, but constrained by transducer's valid span
	    if (transducer.valid_span(first, last, lattice.shortest_distance(first, last))) {
	      // advance dots....
	      
	      // first, extend active items...
	      active_set_type& cell = actives[table](first, last);
	      for (size_t middle = first + 1; middle < last; ++ middle) {
		const active_set_type&  active_arcs  = actives[table](first, middle);
		const passive_set_type& passive_arcs = passives(middle, last);
		
		extend_actives(transducer, active_arcs, passive_arcs, cell);
	      }
	      
	      // then, advance by terminal(s) at lattice[last - 1];
	      const active_set_type&  active_arcs  = actives[table](first, last - 1);
	      const lattice_type::arc_set_type& passive_arcs = lattice[last - 1];
	      
	      active_set_type::const_iterator aiter_begin = active_arcs.begin();
	      active_set_type::const_iterator aiter_end = active_arcs.end();
	      
	      if (aiter_begin != aiter_end) {
		lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
		for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) {
		  const symbol_type& terminal = piter->label;
		  const int length = piter->distance;
		  
		  active_set_type& cell = actives[table](first, last - 1 + length);
		  
		  // handling of EPSILON rule...
		  if (terminal == vocab_type::EPSILON) {
		    for (active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
		      cell.push_back(active_type(aiter->node, aiter->tails, aiter->features + piter->features, aiter->attributes));
		  } else {
		    for (active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
		      const transducer_type::id_type node = transducer.next(aiter->node, terminal);
		      if (node == transducer.root()) continue;
		      
		      cell.push_back(active_type(node, aiter->tails, aiter->features + piter->features, aiter->attributes));
		    }
		  }
		}
	      }
	    }
	    
	    // complete active items if possible... The active items may be created from child span due to the
	    // lattice structure...
	    // apply rules on actives at [first, last)
	    
	    active_set_type&  cell         = actives[table](first, last);
	    passive_set_type& passive_arcs = passives(first, last);
	    
	    active_set_type::const_iterator citer_end = cell.end();
	    for (active_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) {
	      const transducer_type::rule_pair_set_type& rules = transducer.rules(citer->node);
	      
	      if (rules.empty()) continue;
	      
	      transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	      for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
		apply_rule(yield_source ? riter->source : riter->target, riter->features + citer->features, riter->attributes + citer->attributes,
			   citer->tails.begin(), citer->tails.end(), node_map, passive_arcs, graph,
			   first, last);
	    }
	  }
	  
	  // handle unary rules...
	  // order is very important...
	  // actually, we need to loop-forever...

	  if (! passives(first, last).empty()) {
	    passive_set_type& passive_arcs = passives(first, last);
	    
	    size_t passive_first = 0;
	    
	    closure.clear();
	    passive_set_type::const_iterator piter_end = passive_arcs.end();
	    for (passive_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter)
	      closure.insert(non_terminals[*piter]);
	    
	    // run 4 iterations... actually, we should loop until convergence which will be impractical.
	    int  closure_loop = 0;
	    for (int level = 0; /**/; ++ level) {
	      const size_t passive_size = passive_arcs.size();
	      const size_t closure_size = closure.size();
	      
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
		  
		  transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
		  for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
		    const rule_ptr_type rule = (yield_source ? riter->source : riter->target);
		    
		    apply_rule(rule, riter->features, riter->attributes,
			       &passive_arcs[p], (&passive_arcs[p]) + 1, node_map, passive_arcs, graph,
			       first, last, level + 1);
		    
		    closure.insert(rule->lhs);
		  }
		}
	      }
	      
	      if (passive_size == passive_arcs.size()) break;
	      
	      passive_first = passive_size;

	      if (closure_size != closure.size())
		closure_loop = 0;
	      else
		++ closure_loop;
	      
	      if (closure_loop == 4) break;
	    }
	  }
	  
	  // extend root with passive items at [first, last)
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	    
	    if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
	    
	    const active_set_type&  active_arcs  = actives[table](first, first);
	    const passive_set_type& passive_arcs = passives(first, last);
	    
	    active_set_type& cell = actives[table](first, last);
	    
	    extend_actives(transducer, active_arcs, passive_arcs, cell);
	  }
	}
      
      // finally, collect all the parsed rules, and proceed to [goal] rule...
      // passive arcs will not be updated!
      
      // we will clear node map so that we will always create new node..
      node_map.clear();
      
      passive_set_type& passive_arcs = passives(0, lattice.size());
      for (size_t p = 0; p != passive_arcs.size(); ++ p)
	if (non_terminals[passive_arcs[p]] == goal)
	  apply_rule(goal_rule, feature_set_type(), attribute_set_type(), &(passive_arcs[p]), (&passive_arcs[p]) + 1, node_map, passive_arcs, graph,
		     0, lattice.size());
      
      // we will sort to remove unreachable nodes......
      graph.topologically_sort();
    }

  private:
    
    template <typename Iterator>
    void apply_rule(const rule_ptr_type& rule,
		    const feature_set_type& features,
		    const attribute_set_type& attributes,
		    Iterator first,
		    Iterator last,
		    node_map_type& node_map,
		    passive_set_type& passives,
		    hypergraph_type& graph,
		    const int lattice_first,
		    const int lattice_last,
		    const int level = 0)
    {
      hypergraph_type::edge_type& edge = graph.add_edge(first, last);
      edge.rule = rule;
      edge.features = features;
      edge.attributes = attributes;
      
      // assign metadata...
      edge.attributes[attr_span_first] = attribute_set_type::int_type(lattice_first);
      edge.attributes[attr_span_last]  = attribute_set_type::int_type(lattice_last);
      
      std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(rule->lhs, level), 0));
      if (result.second) {
	hypergraph_type::node_type& node = graph.add_node();
	
	if (rule->lhs == goal_rule->lhs)
	  graph.goal = node.id;
	else
	  passives.push_back(node.id);
	
	if (node.id >= non_terminals.size())
	  non_terminals.resize(node.id + 1);
	non_terminals[node.id] = rule->lhs;
	
	result.first->second = node.id;
      }
      
      graph.connect_edge(edge.id, result.first->second);

#if 0
      std::cerr << "new rule: " << *(edge.rule)
		<< " head: " << edge.head
		<< ' ';
      std::copy(edge.tails.begin(), edge.tails.end(), std::ostream_iterator<int>(std::cerr, " "));
      std::cerr << std::endl;
#endif
      
    }
    
    bool extend_actives(const transducer_type& transducer,
			const active_set_type& actives, 
			const passive_set_type& passives,
			active_set_type& cell)
    {
      active_set_type::const_iterator aiter_begin = actives.begin();
      active_set_type::const_iterator aiter_end = actives.end();
      
      passive_set_type::const_iterator piter_begin = passives.begin();
      passive_set_type::const_iterator piter_end = passives.end();

      bool found = false;
      
      if (piter_begin != piter_end)
	for (active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  if (transducer.has_next(aiter->node))
	    for (passive_set_type::const_iterator piter = piter_begin; piter != piter_end; ++ piter) {
	      const symbol_type& non_terminal = non_terminals[*piter];
	      
	      const transducer_type::id_type node = transducer.next(aiter->node, non_terminal);
	      if (node == transducer.root()) continue;
	      
	      hypergraph_type::edge_type::node_set_type tails(aiter->tails.size() + 1, *piter);
	      std::copy(aiter->tails.begin(), aiter->tails.end(), tails.begin());
	      
	      cell.push_back(active_type(node, tails, aiter->features, aiter->attributes));
	      
	      found = true;
	    }
      
      return found;
    }


    
  private:
    const symbol_type goal;
    const grammar_type& grammar;
    const bool yield_source;
    const attribute_type attr_span_first;
    const attribute_type attr_span_last;
    
    rule_ptr_type goal_rule;

    active_chart_set_type  actives;
    passive_chart_type     passives;

    node_map_type         node_map;
    closure_type          closure;
    non_terminal_set_type non_terminals;
  };
  
  inline
  void compose_cky(const Symbol& goal, const Grammar& grammar, const Lattice& lattice, HyperGraph& graph, const bool yield_source=false)
  {
    ComposeCKY(goal, grammar, yield_source)(lattice, graph);
  }
};

#endif
