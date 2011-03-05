// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_EARLEY__HPP__
#define __CICADA__COMPOSE_EARLEY__HPP__ 1


#include <vector>
#include <deque>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/sort.hpp>

#include <google/dense_hash_set>
#include <google/dense_hash_map>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/hashmurmur.hpp>

namespace cicada
{
  
  struct ComposeEarley
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;

    typedef Lattice    lattice_type;
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    ComposeEarley(const grammar_type& __grammar, const bool __yield_source=false)
      : grammar(__grammar), yield_source(__yield_source)

    {
      items_unique_active.set_empty_key(0);
      items_unique_passive.set_empty_key(0);

      traversals.set_empty_key(traversal_type());

      terminal_nodes.set_empty_key(transducer_id_type(-1, 0));
      non_terminal_nodes.set_empty_key(0);
      
      goal_nodes.set_empty_key(hypergraph_type::id_type(-1));
    }
    
    //
    // compose source hypergraph with FST grammar in grammar...!
    //
    

    typedef uint32_t id_type;
    typedef google::dense_hash_map<symbol_type, id_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > id_map_type;
    
  
    // we assume that we have only unique path from tail-nodes to head-node...
    struct grammar_node_type
    {
      grammar_node_type(const bool __is_root)
	: edge(hypergraph_type::invalid), is_root(__is_root) { initialize(); }
      grammar_node_type()
	: edge(hypergraph_type::invalid), is_root(false) { initialize(); }
      
      id_map_type terminals;
      id_map_type non_terminals;
      
      hypergraph_type::id_type edge;
      bool is_root;

    private:
      void initialize()
      {
	terminals.set_empty_key(symbol_type());
	non_terminals.set_empty_key(symbol_type());
      }
    };

    typedef utils::chunk_vector<grammar_node_type, 4096 / sizeof(grammar_node_type), std::allocator<grammar_node_type> > grammar_node_set_type;

    typedef std::pair<int, transducer_type::id_type> transducer_id_type;

    struct Item
    {
      typedef Item item_type;

      symbol_type lhs; // lhs of rule
      const grammar_node_type* dot;  // earley grammar position
      
      // transducer's start and end node span
      transducer_id_type first;
      transducer_id_type last;
      
      // transducer's node for terminal
      transducer_id_type terminal;
      
      // backpointers
      const item_type* active;
      const item_type* passive;
      
      // backptr to source hypergraph's edge
      hypergraph_type::id_type edge;
      
      // created by predict...
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& q0)
	: lhs(__lhs), dot(&__dot),
	  first(q0), last(q0), terminal(transducer_id_type(-1, 0)),
	  active(0), passive(0),
	  edge(hypergraph_type::invalid) {}
      
      // created by predict...
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& q0,
	   const item_type& __active)
	: lhs(__lhs), dot(&__dot),
	  first(q0), last(q0), terminal(transducer_id_type(-1, 0)),
	  active(&__active), passive(0),
	  edge(hypergraph_type::invalid) {}

      
      // created by scan... we will always have terminal
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& __first, const transducer_id_type& __last, const transducer_id_type& __terminal,
	   const item_type& __active,
	   const hypergraph_type::id_type& __edge)
	: lhs(__lhs), dot(&__dot),
	  first(__first), last(__last), terminal(__terminal),
	  active(&__active), passive(0),
	  edge(__edge) {}
      
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& __first, const transducer_id_type& __last, const transducer_id_type& __terminal,
	   const item_type& __active)
	: lhs(__lhs), dot(&__dot),
	  first(__first), last(__last), terminal(__terminal),
	  active(&__active), passive(0),
	  edge(hypergraph_type::invalid) {}
      
      // construct by complete
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& __first, const transducer_id_type& __last,
	   const item_type& __active, const item_type& __passive,
	   const hypergraph_type::id_type& __edge)
	: lhs(__lhs), dot(&__dot),
	  first(__first), last(__last), terminal(transducer_id_type(-1, 0)),
	  active(&__active), passive(&__passive),
	  edge(__edge) {}
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& __first, const transducer_id_type& __last,
	   const item_type& __active, const item_type& __passive)
	: lhs(__lhs), dot(&__dot),
	  first(__first), last(__last), terminal(transducer_id_type(-1, 0)),
	  active(&__active), passive(&__passive),
	  edge(hypergraph_type::invalid) {}
      
      // for query...
      Item(const transducer_id_type& __first, const transducer_id_type& __last)
	: first(__first), last(__last) {}
      
    public:
      bool is_passive() const { return edge != hypergraph_type::invalid; }
      bool is_active() const { return edge == hypergraph_type::invalid; }
      bool is_scanned() const { return active && !passive && dot != 0; }
      bool is_predicted() const { return dot->is_root; }
      bool is_completed() const { return active && passive; }
    };
    
    typedef Item item_type;


    typedef utils::chunk_vector<item_type, 1024 * 16 / sizeof(item_type), std::allocator<item_type> > item_set_type;
    typedef std::deque<const item_type*, std::allocator<const item_type*> > agenda_type;
    
    
    
    // traversal structure to keep only unique items...
    struct Traversal
    {
      const item_type* active;
      const item_type* passive;
      int is_active;
      
      Traversal(const item_type* __active, const item_type* __passive, const bool __is_active)
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
    
    // item hash/comparison
    
    struct item_unique_active_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const item_type* x) const
      {
	// passive item do not care about what dot we have consumed...
	return (x ? hasher_type::operator()(x->dot, hasher_type::operator()(x->first, hasher_type::operator()(x->last, x->lhs.id()))) : size_t(0));
      }
    };
    
    struct item_unique_passive_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const item_type* x) const
      {
	// passive item do not care about what dot we have consumed...
	return (x ? hasher_type::operator()(x->first, hasher_type::operator()(x->last, x->lhs.id())) : size_t(0));
      }
    };
    
    struct item_unique_active_equal_type
    {
      bool operator()(const item_type* x, const item_type* y) const
      {
	// passive item do not care dot...
	return ((x == y) || (x && y && x->lhs == y->lhs && x->first == y->first && x->last == y->last && x->dot == y->dot));
      }
    };
    struct item_unique_passive_equal_type
    {
      bool operator()(const item_type* x, const item_type* y) const
      {
	// passive item do not care dot...
	return ((x == y) || (x && y && x->lhs == y->lhs && x->first == y->first && x->last == y->last));
      }
    };

    
    struct item_active_hash_type : public utils::hashmurmur<size_t>
    {
      size_t operator()(const item_type* x) const
      {
	return utils::hashmurmur<size_t>::operator()(x->last);
      }
    };
    
    struct item_active_equal_type
    {
      bool operator()(const item_type* x, const item_type* y) const
      {
	return x->last == y->last;
      }
    };

    struct item_passive_hash_type : public utils::hashmurmur<size_t>
    {
      size_t operator()(const item_type* x) const
      {
	return utils::hashmurmur<size_t>::operator()(x->first);
      }
    };
    
    struct item_passive_equal_type
    {
      bool operator()(const item_type* x, const item_type* y) const
      {
	return x->first == y->first;
      }
    };
    
    typedef google::dense_hash_set<const item_type*, item_unique_active_hash_type, item_unique_active_equal_type >   item_set_unique_active_type;
    typedef google::dense_hash_set<const item_type*, item_unique_passive_hash_type, item_unique_passive_equal_type > item_set_unique_passive_type;

#ifdef HAVE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_multiset<const item_type*, item_active_hash_type, item_active_equal_type,
					 std::allocator<const item_type*> > item_set_active_type;
    typedef std::tr1::unordered_multiset<const item_type*, item_passive_hash_type, item_passive_equal_type,
					 std::allocator<const item_type*> > item_set_passive_type;
#else
    typedef sgi::hash_multiset<const item_type*, item_active_hash_type, item_active_equal_type,
			       std::allocator<const item_type*> > item_set_active_type;
    typedef sgi::hash_multiset<const item_type*, item_passive_hash_type, item_passive_equal_type,
			       std::allocator<const item_type*> > item_set_passive_type;
#endif
    
    // item to traversal graph mappings...
    
    typedef google::dense_hash_map<transducer_id_type, hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<transducer_id_type> > terminal_node_set_type;
    typedef google::dense_hash_map<const item_type*, hypergraph_type::id_type, item_unique_hash_type, item_unique_equal_type > non_terminal_node_set_type;
    typedef google::dense_hash_set<hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<hypergraph_type::id_type> > goal_node_set_type;

    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // grammar-node. we will use zero for root..
      // thus, we have not resize to zero everytime.

      target.clear();
      if (source.goal == hypergraph_type::invalid)
	return;
      
      initialize_grammar(source);

      id_map_type::const_iterator giter = grammar_nodes[0].non_terminals.find(goal_symbol);
      if (giter == grammar_nodes[0].non_terminals.end())
	throw std::runtime_error("no grammar to accept goal?");

      const grammar_node_type& dot_goal = grammar_nodes[giter->second];
      
      // initial items for each transducer defined in grammar...
      insert_item(item_type(goal_symbol, dot_goal, transducer_id_type(-1, 0)));
      
      // forever...
      while (! agenda_exploration.empty() || ! agenda_finishing.empty()) {
	
	// connect all the newly created items at previous step.
	agenda_type::const_iterator aiter_end = agenda_exploration.end();
	for (agenda_type::const_iterator aiter = agenda_exploration.begin(); aiter != aiter_end; ++ aiter)
	  connect_item(*(*aiter), source, target);
	
	agenda_exploration.clear();
	
	if (! agenda_finishing.empty()) {
	  const item_type* item = agenda_finishing.front();
	  agenda_finishing.pop_front();
	  
	  if (item->is_active()) {
	    if (! item->dot->terminals.empty())
	      scan(*item);
	    
	    if (! item->dot->non_terminals.empty()) {
	      complete_active(*item);
	      predict(*item);
	    }
	  } else
	    complete_passive(*item);
	}
      }
      
      if (! goal_nodes.empty()) {
	// add new node for root...
	target.goal = target.add_node().id;
	
	goal_node_set_type::const_iterator niter_end = goal_nodes.end();
	for (goal_node_set_type::const_iterator niter = goal_nodes.begin(); niter != niter_end; ++ niter) {
	  hypergraph_type::edge_type& edge = target.add_edge(&(*niter), &(*niter) + 1);
	  
	  edge.rule = rule_goal;
	  
	  target.connect_edge(edge.id, target.goal);
	}
      }
      
      if (target.goal != hypergraph_type::invalid)
	remove_epsilon(target);
    }

  private:
    
    void scan(const item_type& item)
    {
      const grammar_node_type& dot = *item.dot;
      
      if (item.last.first >= 0) {
	const int table = item.last.first;
	
	const transducer_type::id_type last = item.last.second;
	const transducer_type& transducer = grammar[table];
	
	id_map_type::const_iterator titer_end = dot.terminals.end();
	for (id_map_type::const_iterator titer = dot.terminals.begin(); titer != titer_end; ++ titer) {
	  
	  const transducer_type::id_type last_next = transducer.next(last, titer->first);
	  if (last_next == transducer.root()) continue;
	  
	  const grammar_node_type& dot_next = grammar_nodes[titer->second];
	  
	  const bool has_rule = dot_next.edge != hypergraph_type::invalid;
	  const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	  
	  // test if we reached a leaf...
	  if (transducer.has_next(last_next)) {
	    if (has_rule)
	      insert_item(item_type(item.lhs, dot_next, item.first, std::make_pair(table, last_next), transducer_id_type(-1, 0), item, dot_next.edge));
	    if (has_next)
	      insert_item(item_type(item.lhs, dot_next, item.first, std::make_pair(table, last_next), transducer_id_type(-1, 0), item));
	  }
	  
	  // test if we have a phrase... we will back to root()
	  if (! transducer.rules(last_next).empty()) {
	    if (has_rule)
	      insert_item(item_type(item.lhs, dot_next, item.first, transducer_id_type(-1, 0), std::make_pair(table, last_next), item, dot_next.edge));
	    if (has_next)
	      insert_item(item_type(item.lhs, dot_next, item.first, transducer_id_type(-1, 0), std::make_pair(table, last_next), item));
	  }
	}
      } else {
	for (size_t table = 0; table != grammar.size(); ++ table) {
	  const transducer_type& transducer = grammar[table];
	  const transducer_type::id_type last = transducer.root();
	  
	  id_map_type::const_iterator titer_end = dot.terminals.end();
	  for (id_map_type::const_iterator titer = dot.terminals.begin(); titer != titer_end; ++ titer) {
	    
	    const transducer_type::id_type last_next = transducer.next(last, titer->first);
	    if (last_next == transducer.root()) continue;
	    
	    const grammar_node_type& dot_next = grammar_nodes[titer->second];
	    
	    const bool has_rule = dot_next.edge != hypergraph_type::invalid;
	    const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	    
	    // test if we reached a leaf...
	    if (transducer.has_next(last_next)) {
	      if (has_rule)
		insert_item(item_type(item.lhs, dot_next, item.first, std::make_pair(table, last_next), transducer_id_type(-1, 0), item, dot_next.edge));
	      if (has_next)
		insert_item(item_type(item.lhs, dot_next, item.first, std::make_pair(table, last_next), transducer_id_type(-1, 0), item));
	    }
	    
	    // test if we have a phrase... we will back to root()
	    if (! transducer.rules(last_next).empty()) {
	      if (has_rule)
		insert_item(item_type(item.lhs, dot_next, item.first, transducer_id_type(-1, 0), std::make_pair(table, last_next), item, dot_next.edge));
	      if (has_next)
		insert_item(item_type(item.lhs, dot_next, item.first, transducer_id_type(-1, 0), std::make_pair(table, last_next), item));
	    }
	  }
	}
      }
    }
    
    void predict(const item_type& item)
    {
      const grammar_node_type& dot = *item.dot;
      
      id_map_type::const_iterator niter_end = dot.non_terminals.end();
      for (id_map_type::const_iterator niter = dot.non_terminals.begin(); niter != niter_end; ++ niter) {
	
	// we will predict this non-terminal
	const symbol_type& lhs = niter->first;
	
	// lookup lhs...
	id_map_type::const_iterator piter = grammar_nodes.front().non_terminals.find(lhs);
	if (piter == grammar_nodes.front().non_terminals.end())
	  continue;
	
	const grammar_node_type& dot_next = grammar_nodes[piter->second];
	
	insert_item(item_type(lhs, dot_next, item.last, item));
      }
    }

    // comlete passive with actives
    void complete_passive(const item_type& passive)
    {
      // we will try find actives whose last match with passive's first
      const item_type query(passive.first, passive.first);
      
      std::pair<item_set_active_type::const_iterator, item_set_active_type::const_iterator> result = items_active.equal_range(&query);
      for (item_set_active_type::const_iterator aiter = result.first; aiter != result.second; ++ aiter) {
	const item_type& active = *(*aiter);
	
	const grammar_node_type& dot = *(active.dot);
	id_map_type::const_iterator niter = dot.non_terminals.find(passive.lhs);
	if (niter == dot.non_terminals.end()) continue;
	
	const grammar_node_type& dot_next = grammar_nodes[niter->second];
	
	const bool has_rule = dot_next.edge != hypergraph_type::invalid;
	const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	
	if (has_rule)
	  insert_item(item_type(active.lhs, dot_next, active.first, passive.last, active, passive, dot_next.edge));
	if (has_next)
	  insert_item(item_type(active.lhs, dot_next, active.first, passive.last, active, passive));
      }
    }
    
    // find passives that can extend the active
    void complete_active(const item_type& active)
    {
      const grammar_node_type& dot = *(active.dot);
      
      // find passives whose first match with active's last
      const item_type query(active.last, active.last);
      
      std::pair<item_set_passive_type::const_iterator, item_set_passive_type::const_iterator> result = items_passive.equal_range(&query);
      for (item_set_passive_type::const_iterator piter = result.first; piter != result.second; ++ piter) {
	const item_type& passive = *(*piter);
	
	id_map_type::const_iterator niter = dot.non_terminals.find(passive.lhs);
	if (niter == dot.non_terminals.end()) continue;

	const grammar_node_type& dot_next = grammar_nodes[niter->second];
	
	const bool has_rule = dot_next.edge != hypergraph_type::invalid;
	const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	
	if (has_rule)
	  insert_item(item_type(active.lhs, dot_next, active.first, passive.last, active, passive, dot_next.edge));
	if (has_next)
	  insert_item(item_type(active.lhs, dot_next, active.first, passive.last, active, passive));
      }
    }

    
    void insert_item(const item_type& item)
    {
      if (item.passive && item.active)
	if (! traversals.insert(traversal_type(item.active, item.passive, item.is_active())).second)
	  return;
      
      items.push_back(item);
      
      agenda_exploration.push_back(&items.back());
    }
    
    void connect_item(const item_type& item, const hypergraph_type& source, hypergraph_type& target)
    {
      if (items.is_passive()) {
	if (items_unique_passive.insert(&item).second) {
	  items_passive.insert(&item);
	  agenda_finishing.push_back(&item);
	}
      } else {
	if (items_unique_active.insert(&item).second) {
	  items_active.insert(&item);
	  agenda_finishing.push_back(&item);
	}
      }
      
      // add into hypergraph...

      const hypergraph_type::node_type* terminal = 0;
      
      if (item.terminal.first >= 0) {
	const grammar_type::rule_pair_set_type& rules = grammar[item.terminal.first].rules(item.terminal.second);
	
	if (! rules.empty()) {
	  
	  terminal_node_set_type::iterator niter = terminal_nodes.find(item.terminal);
	  if (niter == terminal_nodes.end()) {
	    
	    hypergraph_type::node_type& node = target.add_node();
	    
	    grammar_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	    for (grammar_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	      
	      hypergraph_type::edge_type& edge = target.add_edge();
	      edge.rule = (yield_source ? riter->source : riter->target);
	      edge.features   = riter->features;
	      edge.attributes = riter->attributes;
	      
	      target.connect_edge(edge.id, node.id);
	    }
	    
	    niter = terminal_nodes.insert(std::make_pair(item.terminal, node.id)).first;
	  }
	  
	  terminal = &target.nodes[niter->second];
	}
      }
      
      // now, get node...
      non_terminal_node_set_type::iterator niter = non_terminal_nodes.find(&item);
      if (niter == non_terminal_nodes.end())
	niter = non_terminal_nodes.insert(std::make_pair(&item, target.add_node().id)).first;
      
      hypergraph_type::node_type& node_head = target.nodes[niter->second];
      
      if (item.lhs == goal_symbol
	  && item.is_passive()
	  && item.first == transducer_id_type(-1, 0)
	  && item.last == transducer_id_type(-1, 0))
	goal_nodes.insert(node_head.id);
      
      std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails;
      
      if (item.is_predicted()) {
	
      } else if (item.is_scanned()) {
	// we will have one node... + terminal node if terminal exists!

	non_terminal_node_set_type::iterator niter = non_terminal_nodes.find(item.active);
	if (niter == non_terminal_nodes.end())
	  throw std::runtime_error("error during scanning?");
	
	tails.push_back(niter->second);
	if (terminal)
	  tails.push_back(terminal->id);
	
      } else if (item.is_completed()) {
	// we will have two nodes... active and passive
	
	non_terminal_node_set_type::const_iterator niter_active = non_terminal_nodes.find(item.active);
	if (niter_active == non_terminal_nodes.end())
	  throw std::runtime_error("error during completion for active?");

	non_terminal_node_set_type::const_iterator niter_passive = non_terminal_nodes.find(item.passive);
	if (niter_passive == non_terminal_nodes.end())
	  throw std::runtime_error("error during completion for passive?");
	
	tails.push_back(niter_active->second);
	tails.push_back(niter_passive->second);
      } else
	throw std::runtime_error("where this item comes from?");
      
      hypergraph_type::edge_type& graph_edge = target.add_edge(tails.begin(), tails.end());
      
      // TODO: add rule
      switch (tails.size()) {
      case 0: graph_edge.rule = rule_epsilon; break;
      case 1: graph_edge.rule = rule_x1;      break;
      case 2: graph_edge.rule = rule_x1_x2;   break;
      }

      if (item.edge != hypergraph_type::invalid) {
	graph_edge.features   = source.edges[item.edge].features;
	graph_edge.attributes = source.edges[item.edge].attributes;
      }
      
      target.connect_edge(graph_edge.id, node_head.id);
    }


    void initialize_grammar(const hypergraph_type& source)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
      
      items.clear();
      
      traversals.clear();
      terminal_nodes.clear();
      non_terminal_nodes.clear();

      agenda_finishing.clear();
      agenda_exploration.clear();

      items_unique_active.clear();
      items_unique_passive.clear();
      items_active.clear();
      items_passive.clear();

      goal_nodes.clear();
      
      if (! rule_epsilon)
	rule_epsilon = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(1, vocab_type::EPSILON)));
      
      if (! rule_goal)
	rule_goal = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, vocab_type::X1)));
      
      if (! rule_x1)
	rule_x1 = rule_type::create(rule_type(vocab_type::X, rule_type::symbol_set_type(1, vocab_type::X1)));
      
      if (! rule_x1_x2) {
	std::vector<symbol_type, std::allocator<symbol_type> > sequence(2);
	sequence.front() = vocab_type::X1;
	sequence.back() = vocab_type::X2;
	
	rule_x1_x2 =  rule_type::create(rule_type(vocab_type::X, sequence.begin(), sequence.end()));
      }

      
      // assigne pseudo non-terminals
      non_terminal_set_type non_terminals(source.nodes.size());
      for (size_type id = 0; id < source.nodes.size(); ++ id) {
	non_terminals[id] = std::string("[NODE_") + boost::lexical_cast<std::string>(id) + ']';
	//non_terminals[id] = source.edges[source.nodes[id].edges.front()].rule->lhs.non_terminal();
      }

      // assign goal-symbol!
      goal_symbol = non_terminals[source.goal];
      
      
      
      grammar_nodes.clear();
      grammar_nodes.resize(1);
      
      for (size_type id = 0; id < source.edges.size(); ++ id) {
	const hypergraph_type::edge_type& edge = source.edges[id];
	const symbol_type& non_terminal = non_terminals[edge.head];
	
	id_map_type::iterator niter = grammar_nodes[0].non_terminals.find(non_terminal);
	if (niter == grammar_nodes[0].non_terminals.end()) {
	  const id_type id = grammar_nodes.size();
	  grammar_nodes.push_back(grammar_node_type(true));
	  niter = grammar_nodes[0].non_terminals.insert(std::make_pair(non_terminal, id)).first;
	}
	
	int tail_pos = 0;
	rule_type::symbol_set_type::const_iterator siter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator siter = edge.rule->rhs.begin(); siter != siter_end; ++ siter) {
	  if (siter->is_non_terminal()) {
	    const symbol_type& non_terminal = non_terminals[edge.tails[tail_pos]];
	    
	    const id_type grammar_node = niter->second;
	    
	    niter = grammar_nodes[grammar_node].non_terminals.find(non_terminal);
	    if (niter == grammar_nodes[grammar_node].non_terminals.end()) {
	      const id_type id = grammar_nodes.size();
	      grammar_nodes.push_back(grammar_node_type());
	      niter = grammar_nodes[grammar_node].non_terminals.insert(std::make_pair(non_terminal, id)).first;
	    }
	    
	    ++ tail_pos;
	  } else {
	    const symbol_type& terminal = *siter;
	    
	    const id_type grammar_node = niter->second;
	    
	    niter = grammar_nodes[grammar_node].terminals.find(terminal);
	    if (niter == grammar_nodes[grammar_node].terminals.end()) {
	      const id_type id = grammar_nodes.size();
	      grammar_nodes.push_back(grammar_node_type());
	      niter = grammar_nodes[grammar_node].terminals.insert(std::make_pair(terminal, id)).first;
	    }
	  }
	}
	
	grammar_nodes[niter->second].edge = edge.id;
      }
    }


    typedef std::vector<bool, std::allocator<bool> > removed_type;

    struct filter_epsilon
    {
      const removed_type& removed;
      
      filter_epsilon(const removed_type& __removed) : removed(__removed) {}
      
      template <typename Edge>
      bool operator()(const Edge& edge) const
      {
	return removed[edge.id];
      }
    };
    

    void remove_epsilon(hypergraph_type& graph)
    {
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;
      typedef std::vector<edge_set_type, std::allocator<edge_set_type> > node_map_type;

      hypergraph_type graph_removed;

      removed_type removed(graph.edges.size(), false);
      node_map_type out_edges(graph.nodes.size());
      
      bool new_epsilon = false;
      
      do {
	removed.clear();
	removed.resize(graph.edges.size(), false);

	// compute out-edges...
	{
	  out_edges.clear();
	  out_edges.resize(graph.nodes.size());

	  hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
	  for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	    const hypergraph_type::edge_type& edge = *eiter;
	    
	    hypergraph_type::edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	    for (hypergraph_type::edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	      out_edges[*niter].push_back(edge.id);
	  }
	}
	
	// then, iterate again... and mark removed for epsilon rules...
	hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
	for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter)
	  if (eiter->tails.empty() && eiter->rule == rule_epsilon) {
	    const hypergraph_type::edge_type& edge = *eiter;
	    
	    removed[edge.id] = true;
	    
	    if (! edge.features.empty()) {
	      if (graph.nodes[edge.head].edges.size() != 1) {
		// cannot propagate features...
		// or, do we evenly distribute this...?
	      } else {
		edge_set_type::const_iterator eiter_end = out_edges[edge.head].end();
		for (edge_set_type::const_iterator eiter = out_edges[edge.head].begin(); eiter != eiter_end; ++ eiter)
		  graph.edges[*eiter].features += edge.features;
	      }
	    }
	  }
	
	// removed epsilon edges...
	graph_removed.clear();
	topologically_sort(graph, graph_removed, filter_epsilon(removed), false);
	graph.swap(graph_removed);
	
	// re-compute out-edges
	{
	  out_edges.clear();
	  out_edges.resize(graph.nodes.size());
	  
	  hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
	  for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	    const hypergraph_type::edge_type& edge = *eiter;
	    
	    hypergraph_type::edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	    for (hypergraph_type::edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	      out_edges[*niter].push_back(edge.id);
	  }
	}

	new_epsilon = false;
	
	// fix no-edge nodes...
	hypergraph_type::node_set_type::iterator niter_end = graph.nodes.end();
	for (hypergraph_type::node_set_type::iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) 
	  if (niter->edges.empty()) {
	    const hypergraph_type::node_type& node = *niter;
	    
	    edge_set_type::const_iterator eiter_end = out_edges[node.id].end();
	    for (edge_set_type::const_iterator eiter = out_edges[node.id].begin(); eiter != eiter_end; ++ eiter) {
	      hypergraph_type::edge_type& edge = graph.edges[*eiter];
	      
	      if (edge.tails.size() == 2) {
		edge.rule = rule_x1;
		
		hypergraph_type::edge_type::node_set_type tails(1, edge.tails.front() == node.id ? edge.tails.back() : edge.tails.front());
		
		edge.tails = tails;
		
	      } else {
		edge.rule = rule_epsilon;
		edge.tails = hypergraph_type::edge_type::node_set_type();
		
		new_epsilon = true;
	      }
	    }
	  }
	
      } while (new_epsilon);
      
      // final sort...
      graph.topologically_sort();
    }
    
    
  private:  
    const grammar_type& grammar;
    const bool yield_source;
    
    symbol_type           goal_symbol;
    grammar_node_set_type grammar_nodes;

    rule_ptr_type rule_goal;
    rule_ptr_type rule_epsilon;
    rule_ptr_type rule_x1;
    rule_ptr_type rule_x1_x2;

    item_set_type items;

    traversal_set_type traversals;
    terminal_node_set_type     terminal_nodes;
    non_terminal_node_set_type non_terminal_nodes;

    agenda_type agenda_finishing;
    agenda_type agenda_exploration;
    
    item_set_unique_active_type  items_unique_active;
    item_set_unique_passive_type items_unique_passive;
    item_set_active_type  items_active;
    item_set_passive_type items_passive;

    goal_node_set_type goal_nodes;
  };
  
  inline
  void compose_earley(const Grammar& grammar, const HyperGraph& source, HyperGraph& target, const bool yield_source=false)
  {
    ComposeEarley composer(grammar, yield_source);
      
    composer(source, target);
  }
  
};

#endif
