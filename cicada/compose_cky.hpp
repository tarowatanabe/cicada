// -*- mode: c++ -*-

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
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    
    ComposeCKY(const symbol_type& __goal, const grammar_type& __grammar)
      : goal(__goal), grammar(__grammar)
    {
      goal_rule.reset(new rule_type(vocab_type::GOAL,
				    rule_type::symbol_set_type(1, goal.non_terminal(1))));
    }
    
    struct ActiveItem
    {
      ActiveItem(const transducer_type::id_type& __node,
		 const hypergraph_type::edge_type::node_set_type __tails,
		 const feature_set_type& __features)
	: node(__node),
	  tails(__tails),
	  features(__features) {}
      ActiveItem(const transducer_type::id_type& __node,
		 const feature_set_type& __features)
	: node(__node),
	  tails(),
	  features(__features) {}
      ActiveItem(const transducer_type::id_type& __node)
	: node(__node),
	  tails(),
	  features() {}
      
      transducer_type::id_type                  node;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type                          features;
    };
    
    typedef ActiveItem active_type;
    typedef utils::chunk_vector<active_type, 4096 / sizeof(active_type), std::allocator<active_type> > active_set_type;
    typedef utils::chart<active_set_type, std::allocator<active_set_type> > active_chart_type;
    typedef std::vector<active_chart_type, std::allocator<active_chart_type> > active_chart_set_type;

    typedef hypergraph_type::id_type passive_type;
    typedef std::vector<passive_type, std::allocator<passive_type> > passive_set_type;
    typedef utils::chart<passive_set_type, std::allocator<passive_set_type> > passive_chart_type;
    
    struct NodeMap
    {
      typedef google::dense_hash_map<symbol_type, hypergraph_type::id_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > node_map_type;
      
      typedef node_map_type::value_type     value_type;
      
      typedef node_map_type::const_iterator const_iterator;
      typedef node_map_type::iterator       iterator;
      
      NodeMap() : node_map() { node_map.set_empty_key(symbol_type()); }

      inline       iterator find(const symbol_type& key)       { return node_map.find(key); }
      inline const_iterator find(const symbol_type& key) const { return node_map.find(key); }
      
      inline       iterator begin()       { return node_map.begin(); }
      inline const_iterator begin() const { return node_map.begin(); }

      inline       iterator end()       { return node_map.end(); }
      inline const_iterator end() const { return node_map.end(); }
      
      std::pair<iterator, bool> insert(const value_type& x) { return node_map.insert(x); }
      
      node_map_type node_map;
    };
    typedef NodeMap node_map_type;

    typedef utils::chart<node_map_type, std::allocator<node_map_type> > node_map_chart_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;

    typedef std::pair<int, feature_set_type> pos_feature_type;
    typedef std::deque<pos_feature_type, std::allocator<pos_feature_type> > pos_feature_set_type;
    typedef std::vector<pos_feature_set_type, std::allocator<pos_feature_set_type> > epsilon_set_type;
    
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      if (lattice.empty())
	return;
      
      // initialize internal structure...
      epsilons.clear();
      actives.clear();
      passives.clear();
      nodes.clear();
      non_terminals.clear();
      
      epsilons.resize(lattice.size() + 1);
      actives.resize(grammar.size(), active_chart_type(lattice.size() + 1));
      passives.resize(lattice.size() + 1);
      nodes.resize(lattice.size() + 1);
      
      // first, insert default positions...
      for (size_t i = 0; i != epsilons.size(); ++ i)
	epsilons[i].push_back(pos_feature_type(i, feature_set_type()));

      // second, construct closure for epsilons...
      for (int first = lattice.size() - 1; first >= 0; -- first) {
	const lattice_type::arc_set_type& arcs = lattice[first];
	
	lattice_type::arc_set_type::const_iterator aiter_end = arcs.end();
	for (lattice_type::arc_set_type::const_iterator aiter = arcs.begin(); aiter != aiter_end; ++ aiter) 
	  if (aiter->label == vocab_type::EPSILON) {
	    const int last = first + aiter->distance;
	    
	    epsilons[first].push_back(pos_feature_type(last, aiter->features));
	    
	    pos_feature_set_type::const_iterator eiter_end = epsilons[last].end();
	    for (pos_feature_set_type::const_iterator eiter = epsilons[last].begin() + 1; eiter != eiter_end; ++ eiter)
	      epsilons[first].push_back(pos_feature_type(eiter->first, eiter->second + aiter->features));
	  }
      }
      
      
      // initialize active chart
      for (size_t table = 0; table != grammar.size(); ++ table) {
	const transducer_type::id_type root = grammar[table].root();
	
	// pos 0...
	pos_feature_set_type::const_iterator eiter_end = epsilons.front().end();
	for (pos_feature_set_type::const_iterator eiter = epsilons.front().begin(); eiter != eiter_end; ++ eiter)
	  if (grammar[table].valid_span(0, eiter->first, 0))
	    actives[table](0, eiter->first).push_back(active_type(root, eiter->second));
	
	for (size_t pos = 1; pos != lattice.size(); ++ pos)
	  if (grammar[table].valid_span(pos, pos, 0))
	    actives[table](pos, pos).push_back(active_type(root));
      }
      
      for (size_t length = 1; length <= lattice.size(); ++ length)
	for (size_t first = 0; first + length <= lattice.size(); ++ first) {
	  const size_t last = first + length;
	  
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	    
	    if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
	    
	    // advance dots....
	    
	    // first, extend active items...
	    active_set_type& cell = actives[table](first, last);
	    
	    for (size_t middle = first + 1; middle < last; ++ middle) {
	      const active_set_type&  active_arcs  = actives[table](first, middle);
	      const passive_set_type& passive_arcs = passives(middle, last);
	      
	      extend_actives(transducer, active_arcs, passive_arcs, cell, true);
	    }
	    
	    // then, advance by terminal(s) at lattice[last - 1];
	    {
	      const active_set_type&  active_arcs  = actives[table](first, last - 1);
	      const lattice_type::arc_set_type& passive_arcs = lattice[last - 1];

	      active_set_type::const_iterator aiter_begin = active_arcs.begin();
	      active_set_type::const_iterator aiter_end = active_arcs.end();
	      
	      lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
	      for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) 
		if (piter->label != vocab_type::EPSILON) {
		  const symbol_type& terminal = piter->label;
		  const int length = piter->distance;
		  
		  const pos_feature_set_type& eps = epsilons[last - 1 + length];
		  
		  for (active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
		    const transducer_type::id_type node = transducer.next(aiter->node, terminal);
		    if (node == transducer.root()) continue;
		    
		    pos_feature_set_type::const_iterator eiter_end = eps.end();
		    for (pos_feature_set_type::const_iterator eiter = eps.begin(); eiter != eiter_end; ++ eiter) {
		      
		      actives[table](first, eiter->first).push_back(active_type(node, aiter->tails, aiter->features + piter->features + eiter->second));
		    }
		  }
		}
	    }
	    
	    // apply rules on actives at [first, last)
	    node_map_type&    node_map     = nodes(first, last);
	    passive_set_type& passive_arcs = passives(first, last);
	    
	    active_set_type::const_iterator citer_end = cell.end();
	    for (active_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) {
	      const transducer_type::rule_pair_set_type& rules = transducer.rules(citer->node);
	      
	      if (rules.empty()) continue;
	      
	      transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	      for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
		apply_rule(riter->target, riter->features + citer->features, citer->tails.begin(), citer->tails.end(), node_map, passive_arcs, graph,
			   first, last);
	    }
	  }
	  
	  // handle unary rules...
	  // order is very important...
	  
	  passive_set_type& passive_arcs = passives(first, last);
	  node_map_type& node_map        = nodes(first, last);
	  
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	    
	    if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
	    
	    for (size_t p = 0; p != passive_arcs.size(); ++ p) {
	      const symbol_type& non_terminal = non_terminals[passive_arcs[p]];
	      
	      const transducer_type::id_type node = transducer.next(transducer.root(), non_terminal);
	      if (node == transducer.root()) continue;
	      
	      const transducer_type::rule_pair_set_type& rules = transducer.rules(node);
	      
	      if (rules.empty()) continue;
	      
	      // passive_arcs "MAY" be modified!
	      
	      transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	      for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
		apply_rule(riter->target, riter->features, &passive_arcs[p], (&passive_arcs[p]) + 1, node_map, passive_arcs, graph,
			   first, last);
	    }
	  }
	  
	  // extend applied unary rules...
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	    
	    if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
	    
	    const active_set_type&  active_arcs  = actives[table](first, first);
	    const passive_set_type& passive_arcs = passives(first, last);
	    
	    active_set_type& cell = actives[table](first, last);
	    
	    extend_actives(transducer, active_arcs, passive_arcs, cell, false);
	  }
	}
      
      // finally, collect all the parsed rules, and proceed to [goal] rule...
      // passive arcs will not be updated!
      node_map_type&    node_map     = nodes(0, lattice.size());
      passive_set_type& passive_arcs = passives(0, lattice.size());
      for (size_t p = 0; p != passive_arcs.size(); ++ p)
	if (non_terminals[passive_arcs[p]] == goal)
	  apply_rule(goal_rule, feature_set_type(), &(passive_arcs[p]), (&passive_arcs[p]) + 1, node_map, passive_arcs, graph,
		     0, lattice.size());
      
      // we will sort to remove unreachable nodes......
      graph.topologically_sort();
    }

  private:
        
    template <typename Iterator>
    void apply_rule(const rule_ptr_type& rule,
		    const feature_set_type& features,
		    Iterator first,
		    Iterator last,
		    node_map_type& node_map,
		    passive_set_type& passives,
		    hypergraph_type& graph,
		    const int lattice_first,
		    const int lattice_last)
    {
      hypergraph_type::edge_type& edge = graph.add_edge(first, last);
      edge.rule = rule;
      
      if (! features.empty())
	edge.features += features;
      
      // assign metadata...
      edge.first    = lattice_first;
      edge.last     = lattice_last;
      
      node_map_type::iterator niter = node_map.find(rule->lhs);
      if (niter == node_map.end()) {
	hypergraph_type::node_type& node = graph.add_node();
	
	if (rule->lhs == goal_rule->lhs)
	  graph.goal = node.id;
	else
	  passives.push_back(node.id);
	
	if (node.id >= non_terminals.size())
	  non_terminals.resize(node.id + 1);
	non_terminals[node.id] = rule->lhs;
	
	niter = node_map.insert(std::make_pair(rule->lhs, node.id)).first;
      }
      
      graph.connect_edge(edge.id, niter->second);

#if 0
      std::cerr << "new rule: " << *(edge.rule)
		<< " head: " << edge.head
		<< ' ';
      std::copy(edge.tails.begin(), edge.tails.end(), std::ostream_iterator<int>(std::cerr, " "));
      std::cerr << std::endl;
#endif
      
    }
    
    void extend_actives(const transducer_type& transducer,
			const active_set_type& actives, 
			const passive_set_type& passives,
			active_set_type& cell,
			const bool check_empty)
    {
      active_set_type::const_iterator aiter_begin = actives.begin();
      active_set_type::const_iterator aiter_end = actives.end();
	      
      passive_set_type::const_iterator piter_begin = passives.begin();
      passive_set_type::const_iterator piter_end = passives.end();
      
      for (active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
	
	if (check_empty && aiter->node == transducer.root()) continue;
	
	for (passive_set_type::const_iterator piter = piter_begin; piter != piter_end; ++ piter) {
	  const symbol_type& non_terminal = non_terminals[*piter];
	  
	  const transducer_type::id_type node = transducer.next(aiter->node, non_terminal);
	  if (node == transducer.root()) continue;
	  
	  hypergraph_type::edge_type::node_set_type tails(aiter->tails.size() + 1);
	  std::copy(aiter->tails.begin(), aiter->tails.end(), tails.begin());
	  tails.back() = *piter;
	  
	  cell.push_back(active_type(node, tails, aiter->features));
	}
      }
    }


    
  private:
    const symbol_type goal;
    const grammar_type& grammar;
    
    rule_ptr_type goal_rule;

    epsilon_set_type       epsilons;
    active_chart_set_type  actives;
    passive_chart_type     passives;
    node_map_chart_type    nodes;
    
    non_terminal_set_type non_terminals;
  };
  
  inline
  void compose_cky(const Symbol& goal, const Grammar& grammar, const Lattice& lattice, HyperGraph& graph)
  {
    ComposeCKY(goal, grammar)(lattice, graph);
  }
};

#endif
