// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
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
#include <cicada/sort_topologically.hpp>
#include <cicada/remove_epsilon.hpp>

#include <google/dense_hash_set>
#include <google/dense_hash_map>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/simple_vector.hpp>
#include <utils/small_vector.hpp>
#include <utils/lexical_cast.hpp>

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
      typedef utils::small_vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;
      
      grammar_node_type(const bool __is_root)
	: edges(), is_root(__is_root) { initialize(); }
      grammar_node_type()
	: edges(), is_root(false) { initialize(); }
      
      id_map_type terminals;
      id_map_type non_terminals;
      
      edge_set_type edges;
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
    
    struct Edge
    {
      typedef Edge edge_type;

      symbol_type lhs; // lhs of rule
      const grammar_node_type* dot;  // earley grammar position
      
      // transducer's start and end node span
      transducer_id_type first;
      transducer_id_type last;
      
      // transducer's node for terminal
      transducer_id_type terminal;
      
      // backpointers
      const edge_type* active;
      const edge_type* passive;
      
      // backptr to source hypergraph's edge
      const grammar_node_type* edges;
      
      // created by predict...
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& q0)
	: lhs(__lhs), dot(&__dot),
	  first(q0), last(q0), terminal(transducer_id_type(-1, 0)),
	  active(0), passive(0),
	  edges(0) {}
      
      // created by predict...
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& q0,
	   const edge_type& __active)
	: lhs(__lhs), dot(&__dot),
	  first(q0), last(q0), terminal(transducer_id_type(-1, 0)),
	  active(&__active), passive(0),
	  edges(0) {}

      
      // created by scan... we will always have terminal
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& __first, const transducer_id_type& __last, const transducer_id_type& __terminal,
	   const edge_type& __active,
	   const grammar_node_type& __edges)
	: lhs(__lhs), dot(&__dot),
	  first(__first), last(__last), terminal(__terminal),
	  active(&__active), passive(0),
	  edges(&__edges) {}
      
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& __first, const transducer_id_type& __last, const transducer_id_type& __terminal,
	   const edge_type& __active)
	: lhs(__lhs), dot(&__dot),
	  first(__first), last(__last), terminal(__terminal),
	  active(&__active), passive(0),
	  edges(0) {}
      
      // construct by complete
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& __first, const transducer_id_type& __last,
	   const edge_type& __active, const edge_type& __passive,
	   const grammar_node_type& __edges)
	: lhs(__lhs), dot(&__dot),
	  first(__first), last(__last), terminal(transducer_id_type(-1, 0)),
	  active(&__active), passive(&__passive),
	  edges(&__edges) {}
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const transducer_id_type& __first, const transducer_id_type& __last,
	   const edge_type& __active, const edge_type& __passive)
	: lhs(__lhs), dot(&__dot),
	  first(__first), last(__last), terminal(transducer_id_type(-1, 0)),
	  active(&__active), passive(&__passive),
	  edges(0) {}
      
      // for query...
      Edge(const transducer_id_type& __first, const transducer_id_type& __last)
	: first(__first), last(__last) {}
      
    public:
      bool is_passive() const { return edges; }
      bool is_active() const { return ! edges; }
      bool is_scanned() const { return active && ! passive && dot; }
      bool is_predicted() const { return dot->is_root; }
      bool is_completed() const { return active && passive; }
    };
    
    typedef Edge edge_type;


    typedef utils::chunk_vector<edge_type, 1024 * 16 / sizeof(edge_type), std::allocator<edge_type> > edge_set_type;
    typedef std::deque<const edge_type*, std::allocator<const edge_type*> > agenda_type;
    
    // traversal structure to keep only unique edges...
    struct Traversal
    {
      const edge_type* active;
      const edge_type* passive;
      size_type is_active;
      
      Traversal(const edge_type* __active, const edge_type* __passive, const bool __is_active)
	: active(__active), passive(__passive), is_active(__is_active) {}
      
      Traversal()
	: active(0), passive(0), is_active(0) {}
    };
    typedef Traversal traversal_type;
    
    struct traversal_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const traversal_type& x) const
      {
	return hasher_type::operator()(x.active, hasher_type::operator()(x.passive, x.is_active));
      }
    };
    
    struct traversal_equal_type
    {
      bool operator()(const traversal_type& x, const traversal_type& y) const
      {
	return x.passive == y.passive && x.active == y.active && x.is_active == y.is_active;
      }
    };
    
    typedef google::dense_hash_set<traversal_type, traversal_hash_type, traversal_equal_type > traversal_set_type;
    
    // edge hash/comparison
   struct edge_unique_hash_type : public utils::hashmurmur<size_t>
   {
     typedef utils::hashmurmur<size_t> hasher_type;
     
     size_t operator()(const edge_type* x) const
     {
       if (! x)
         return 0;
       else if (x->is_active())
         return hasher_type::operator()(x->dot, hasher_type::operator()(x->first, hasher_type::operator()(x->last, x->lhs.id())));
       else
         return hasher_type::operator()(x->first, hasher_type::operator()(x->last, x->lhs.id()));
     }
   };
    
    struct edge_unique_equal_type
    {
      bool operator()(const edge_type* x, const edge_type* y) const
      {
        // passive edge do not care dot...
	return ((x == y) 
		|| (x && y
		    && x->is_active() == y->is_active()
		    && x->lhs == y->lhs
		    && x->first == y->first
		    && x->last == y->last
		    && (x->is_passive() || x->dot == y->dot)));
      }
    };
    

    
    struct edge_active_hash_type : public utils::hashmurmur<size_t>
    {
      size_t operator()(const edge_type* x) const
      {
	return utils::hashmurmur<size_t>::operator()(x->last);
      }
    };
    
    struct edge_active_equal_type
    {
      bool operator()(const edge_type* x, const edge_type* y) const
      {
	return x->last == y->last;
      }
    };

    struct edge_passive_hash_type : public utils::hashmurmur<size_t>
    {
      size_t operator()(const edge_type* x) const
      {
	return utils::hashmurmur<size_t>::operator()(x->first);
      }
    };
    
    struct edge_passive_equal_type
    {
      bool operator()(const edge_type* x, const edge_type* y) const
      {
	return x->first == y->first;
      }
    };
    

#ifdef HAVE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_multiset<const edge_type*, edge_active_hash_type, edge_active_equal_type,
					 std::allocator<const edge_type*> > edge_set_active_type;
    typedef std::tr1::unordered_multiset<const edge_type*, edge_passive_hash_type, edge_passive_equal_type,
					 std::allocator<const edge_type*> > edge_set_passive_type;
#else
    typedef sgi::hash_multiset<const edge_type*, edge_active_hash_type, edge_active_equal_type,
			       std::allocator<const edge_type*> > edge_set_active_type;
    typedef sgi::hash_multiset<const edge_type*, edge_passive_hash_type, edge_passive_equal_type,
			       std::allocator<const edge_type*> > edge_set_passive_type;
#endif
    
    // edge to traversal graph mappings...
    
    typedef google::dense_hash_map<transducer_id_type, hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<transducer_id_type> > terminal_node_set_type;
    typedef google::dense_hash_map<const edge_type*, hypergraph_type::id_type, edge_unique_hash_type, edge_unique_equal_type > non_terminal_node_set_type;
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
      
      // initial edges for each transducer defined in grammar...
      insert_edge(edge_type(goal_symbol, dot_goal, transducer_id_type(-1, 0)));
      
      // forever...
      while (! agenda_exploration.empty() || ! agenda_finishing.empty()) {
	
	// connect all the newly created edges at previous step.
	agenda_type::const_iterator aiter_end = agenda_exploration.end();
	for (agenda_type::const_iterator aiter = agenda_exploration.begin(); aiter != aiter_end; ++ aiter)
	  connect_edge(*(*aiter), source, target);
	
	agenda_exploration.clear();
	
	if (! agenda_finishing.empty()) {
	  const edge_type* edge = agenda_finishing.front();
	  agenda_finishing.pop_front();
	  
	  if (edge->is_active()) {
	    if (! edge->dot->terminals.empty())
	      scan(*edge);
	    
	    if (! edge->dot->non_terminals.empty()) {
	      complete_active(*edge);
	      predict(*edge);
	    }
	  } else
	    complete_passive(*edge);
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
	cicada::remove_epsilon(target);
    }

  private:
    
    void scan(const edge_type& edge)
    {
      const grammar_node_type& dot = *edge.dot;
      
      if (edge.last.first >= 0) {
	const int table = edge.last.first;
	
	const transducer_type::id_type last = edge.last.second;
	const transducer_type& transducer = grammar[table];
	
	id_map_type::const_iterator titer_end = dot.terminals.end();
	for (id_map_type::const_iterator titer = dot.terminals.begin(); titer != titer_end; ++ titer) {
	  
	  const transducer_type::id_type last_next = transducer.next(last, titer->first);
	  if (last_next == transducer.root()) continue;
	  
	  const grammar_node_type& dot_next = grammar_nodes[titer->second];
	  
	  const bool has_rule = ! dot_next.edges.empty();
	  const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	  
	  // test if we reached a leaf...
	  if (transducer.has_next(last_next)) {
	    if (has_rule)
	      insert_edge(edge_type(edge.lhs, dot_next, edge.first, std::make_pair(table, last_next), transducer_id_type(-1, 0), edge, dot_next));
	    if (has_next)
	      insert_edge(edge_type(edge.lhs, dot_next, edge.first, std::make_pair(table, last_next), transducer_id_type(-1, 0), edge));
	  }
	  
	  // test if we have a phrase... we will back to root()
	  if (! transducer.rules(last_next).empty()) {
	    if (has_rule)
	      insert_edge(edge_type(edge.lhs, dot_next, edge.first, transducer_id_type(-1, 0), std::make_pair(table, last_next), edge, dot_next));
	    if (has_next)
	      insert_edge(edge_type(edge.lhs, dot_next, edge.first, transducer_id_type(-1, 0), std::make_pair(table, last_next), edge));
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
	    
	    const bool has_rule = ! dot_next.edges.empty();
	    const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	    
	    // test if we reached a leaf...
	    if (transducer.has_next(last_next)) {
	      if (has_rule)
		insert_edge(edge_type(edge.lhs, dot_next, edge.first, std::make_pair(table, last_next), transducer_id_type(-1, 0), edge, dot_next));
	      if (has_next)
		insert_edge(edge_type(edge.lhs, dot_next, edge.first, std::make_pair(table, last_next), transducer_id_type(-1, 0), edge));
	    }
	    
	    // test if we have a phrase... we will back to root()
	    if (! transducer.rules(last_next).empty()) {
	      if (has_rule)
		insert_edge(edge_type(edge.lhs, dot_next, edge.first, transducer_id_type(-1, 0), std::make_pair(table, last_next), edge, dot_next));
	      if (has_next)
		insert_edge(edge_type(edge.lhs, dot_next, edge.first, transducer_id_type(-1, 0), std::make_pair(table, last_next), edge));
	    }
	  }
	}
      }
    }
    
    void predict(const edge_type& edge)
    {
      const grammar_node_type& dot = *edge.dot;
      
      id_map_type::const_iterator niter_end = dot.non_terminals.end();
      for (id_map_type::const_iterator niter = dot.non_terminals.begin(); niter != niter_end; ++ niter) {
	
	// we will predict this non-terminal
	const symbol_type& lhs = niter->first;
	
	// lookup lhs...
	id_map_type::const_iterator piter = grammar_nodes.front().non_terminals.find(lhs);
	if (piter == grammar_nodes.front().non_terminals.end())
	  continue;
	
	const grammar_node_type& dot_next = grammar_nodes[piter->second];
	
	insert_edge(edge_type(lhs, dot_next, edge.last, edge));
      }
    }

    // comlete passive with actives
    void complete_passive(const edge_type& passive)
    {
      // we will try find actives whose last match with passive's first
      const edge_type query(passive.first, passive.first);
      
      std::pair<edge_set_active_type::const_iterator, edge_set_active_type::const_iterator> result = edges_active.equal_range(&query);
      for (edge_set_active_type::const_iterator aiter = result.first; aiter != result.second; ++ aiter) {
	const edge_type& active = *(*aiter);
	
	const grammar_node_type& dot = *(active.dot);
	id_map_type::const_iterator niter = dot.non_terminals.find(passive.lhs);
	if (niter == dot.non_terminals.end()) continue;
	
	const grammar_node_type& dot_next = grammar_nodes[niter->second];
	
	const bool has_rule = ! dot_next.edges.empty();
	const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	
	if (has_rule)
	  insert_edge(edge_type(active.lhs, dot_next, active.first, passive.last, active, passive, dot_next));
	if (has_next)
	  insert_edge(edge_type(active.lhs, dot_next, active.first, passive.last, active, passive));
      }
    }
    
    // find passives that can extend the active
    void complete_active(const edge_type& active)
    {
      const grammar_node_type& dot = *(active.dot);
      
      // find passives whose first match with active's last
      const edge_type query(active.last, active.last);
      
      std::pair<edge_set_passive_type::const_iterator, edge_set_passive_type::const_iterator> result = edges_passive.equal_range(&query);
      for (edge_set_passive_type::const_iterator piter = result.first; piter != result.second; ++ piter) {
	const edge_type& passive = *(*piter);
	
	id_map_type::const_iterator niter = dot.non_terminals.find(passive.lhs);
	if (niter == dot.non_terminals.end()) continue;

	const grammar_node_type& dot_next = grammar_nodes[niter->second];
	
	const bool has_rule = ! dot_next.edges.empty();
	const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	
	if (has_rule)
	  insert_edge(edge_type(active.lhs, dot_next, active.first, passive.last, active, passive, dot_next));
	if (has_next)
	  insert_edge(edge_type(active.lhs, dot_next, active.first, passive.last, active, passive));
      }
    }

    
    void insert_edge(const edge_type& edge)
    {
      if (edge.is_completed())
	if (! traversals.insert(traversal_type(edge.active, edge.passive, edge.is_active())).second)
	  return;
      
      edges.push_back(edge);
      
      agenda_exploration.push_back(&edges.back());
    }
    
    void connect_edge(const edge_type& edge, const hypergraph_type& source, hypergraph_type& target)
    {
      // now, get node...
      std::pair<non_terminal_node_set_type::iterator, bool> result = non_terminal_nodes.insert(std::make_pair(&edge, 0));
      if (result.second)
	result.first->second = target.add_node().id;
      
      const hypergraph_type::id_type node_head = result.first->second;
      
      if (result.second) {
	if (edge.is_passive())
	  edges_passive.insert(&edge);
	else
	  edges_active.insert(&edge);
	
	agenda_finishing.push_back(&edge);
      }
      
      // add into hypergraph...

      hypergraph_type::id_type terminal = hypergraph_type::invalid;
      
      if (edge.terminal.first >= 0) {
	const grammar_type::rule_pair_set_type& rules = grammar[edge.terminal.first].rules(edge.terminal.second);
	
	if (! rules.empty()) {
	  
	  terminal_node_set_type::iterator niter = terminal_nodes.find(edge.terminal);
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
	    
	    niter = terminal_nodes.insert(std::make_pair(edge.terminal, node.id)).first;
	  }
	  
	  terminal = niter->second;
	}
      }
            
      if (edge.lhs == goal_symbol
	  && edge.is_passive()
	  && edge.first == transducer_id_type(-1, 0)
	  && edge.last == transducer_id_type(-1, 0))
	goal_nodes.insert(node_head);
      
      std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails;
      
      if (edge.is_predicted()) {
	
      } else if (edge.is_scanned()) {
	// we will have one node... + terminal node if terminal exists!
	
	non_terminal_node_set_type::iterator niter = non_terminal_nodes.find(edge.active);
	if (niter == non_terminal_nodes.end())
	  throw std::runtime_error("error during scanning?");
	
	tails.push_back(niter->second);
	if (terminal != hypergraph_type::invalid)
	  tails.push_back(terminal);
	
      } else if (edge.is_completed()) {
	// we will have two nodes... active and passive
	
	non_terminal_node_set_type::const_iterator niter_active = non_terminal_nodes.find(edge.active);
	if (niter_active == non_terminal_nodes.end()) 
	  throw std::runtime_error("error during completion for active?");
	
	non_terminal_node_set_type::const_iterator niter_passive = non_terminal_nodes.find(edge.passive);
	if (niter_passive == non_terminal_nodes.end())
	  throw std::runtime_error("error during completion for passive?");
	
	tails.push_back(niter_active->second);
	tails.push_back(niter_passive->second);
      } else
	throw std::runtime_error("where this edge comes from?");
      
      if (edge.is_active()) {
	hypergraph_type::edge_type& graph_edge = target.add_edge(tails.begin(), tails.end());
	
	switch (tails.size()) {
	case 0: graph_edge.rule = rule_epsilon; break;
	case 1: graph_edge.rule = rule_x1;      break;
	case 2: graph_edge.rule = rule_x1_x2;   break;
	}
	
	target.connect_edge(graph_edge.id, node_head);
      } else {
	for (size_t i = 0; i != edge.edges->edges.size(); ++ i) {
	  hypergraph_type::edge_type& graph_edge = target.add_edge(tails.begin(), tails.end());
	  
	  switch (tails.size()) {
	  case 0: graph_edge.rule = rule_epsilon; break;
	  case 1: graph_edge.rule = rule_x1;      break;
	  case 2: graph_edge.rule = rule_x1_x2;   break;
	  }
	  
	  graph_edge.features   = source.edges[edge.edges->edges[i]].features;
	  graph_edge.attributes = source.edges[edge.edges->edges[i]].attributes;
	  
	  target.connect_edge(graph_edge.id, node_head);
	}
      }
    }


    void initialize_grammar(const hypergraph_type& source)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
      
      edges.clear();
      
      traversals.clear();
      terminal_nodes.clear();
      non_terminal_nodes.clear();

      agenda_finishing.clear();
      agenda_exploration.clear();

      edges_active.clear();
      edges_passive.clear();

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
	non_terminals[id] = std::string("[NODE_") + utils::lexical_cast<std::string>(id) + ']';
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
	
	int non_terminal_pos = 0;
	rule_type::symbol_set_type::const_iterator siter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator siter = edge.rule->rhs.begin(); siter != siter_end; ++ siter) {
	  if (siter->is_non_terminal()) {
	    const int __non_terminal_index = siter->non_terminal_index();
	    const int pos = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	    ++ non_terminal_pos;
	    
	    const symbol_type& non_terminal = non_terminals[edge.tails[pos]];
	    
	    const id_type grammar_node = niter->second;
	    
	    niter = grammar_nodes[grammar_node].non_terminals.find(non_terminal);
	    if (niter == grammar_nodes[grammar_node].non_terminals.end()) {
	      const id_type id = grammar_nodes.size();
	      grammar_nodes.push_back(grammar_node_type());
	      niter = grammar_nodes[grammar_node].non_terminals.insert(std::make_pair(non_terminal, id)).first;
	    }
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
	
	grammar_nodes[niter->second].edges.push_back(edge.id);
      }
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

    edge_set_type edges;

    traversal_set_type traversals;
    terminal_node_set_type     terminal_nodes;
    non_terminal_node_set_type non_terminal_nodes;

    agenda_type agenda_finishing;
    agenda_type agenda_exploration;
    
    edge_set_active_type  edges_active;
    edge_set_passive_type edges_passive;

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
