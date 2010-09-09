// -*- mode: c++ -*-

#ifndef __CICADA__GENERATE_EARLEY__HPP__
#define __CICADA__GENERATE_EARLEY__HPP__ 1


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
#include <cicada/span_node.hpp>

#include <google/dense_hash_set>
#include <google/dense_hash_map>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/compact_trie.hpp>

namespace cicada
{
  
  struct GenerateEarley
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
    
    GenerateEarley()
    {
      edges_unique.set_empty_key(0);

      traversals.set_empty_key(traversal_type());

      non_terminal_nodes.set_empty_key(0);
    }
    
    //
    // generate source hypergraph with FST grammar in grammar...!
    //
    
    typedef uint32_t id_type;
    typedef google::dense_hash_map<symbol_type, id_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > id_map_type;

    typedef std::pair<int, int> span_type;
    typedef std::vector<span_type, std::allocator<span_type> > span_set_type;
    
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
      feature_set_type features;
      bool is_root;

    private:
      void initialize()
      {
	terminals.set_empty_key(symbol_type());
	non_terminals.set_empty_key(symbol_type());
      }
    };

    typedef utils::chunk_vector<grammar_node_type, 4096 / sizeof(grammar_node_type), std::allocator<grammar_node_type> > grammar_node_set_type;
    
    struct Edge
    {
      typedef Edge edge_type;
      
      symbol_type lhs; // lhs of rule
      const grammar_node_type* dot;  // earley grammar position
      
      // backpointers
      const edge_type* active;
      const edge_type* passive;
      
      // backptr to source hypergraph's edge
      hypergraph_type::id_type edge;
      
      
      // created by predict...
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot)
	: lhs(__lhs), dot(&__dot),
	  active(0), passive(0),
	  edge(hypergraph_type::invalid) {}
      
      // created by predict...
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const edge_type& __active)
	: lhs(__lhs), dot(&__dot),
	  active(&__active), passive(0),
	  edge(hypergraph_type::invalid) {}

      
      // created by scan... we will always have terminal
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const edge_type& __active,
	   const hypergraph_type::id_type& __edge)
	: lhs(__lhs), dot(&__dot),
	  active(&__active), passive(0),
	  edge(__edge) {}
      
      
      // construct by complete
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const edge_type& __active, const edge_type& __passive,
	   const hypergraph_type::id_type& __edge)
	: lhs(__lhs), dot(&__dot),
	  active(&__active), passive(&__passive),
	  edge(__edge) {}
      Edge(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const edge_type& __active, const edge_type& __passive)
	: lhs(__lhs), dot(&__dot),
	  active(&__active), passive(&__passive),
	  edge(hypergraph_type::invalid) {}
      

            
    public:
      bool is_passive() const { return edge != hypergraph_type::invalid; }
      bool is_active() const { return edge == hypergraph_type::invalid; }
      bool is_scanned() const { return active && !passive && dot != 0; }
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
      int is_active;
      
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
    
    // edge hash/comparison
    
    struct edge_unique_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const edge_type* x) const
      {
	// passive edge do not care about what dot we have consumed...
	if (! x)
	  return 0;
	else if (x->is_active())
	  return hasher_type::operator()(x->dot, x->lhs.id());
	else
	  return hasher_type::operator()(x->lhs.id());
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
		    && (x->is_passive() || x->dot == y->dot)));
      }
    };
    
    
    typedef google::dense_hash_set<const edge_type*, edge_unique_hash_type, edge_unique_equal_type > edge_set_unique_type;

    typedef std::vector<const edge_type*, std::allocator<const edge_type*> > edge_set_active_type;
    typedef std::vector<const edge_type*, std::allocator<const edge_type*> > edge_set_passive_type;
    
    // edge to traversal graph mappings...
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
      insert_edge(edge_type(goal_symbol, dot_goal));
      
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
      
      // sort!
      if (target.is_valid())
	target.topologically_sort();
    }

  private:
    
    void scan(const edge_type& edge)
    {
      const grammar_node_type& dot = *edge.dot;

      //
      // scanning implies moving dot on terminals w/o intersecting with transducer...
      //
      
      id_map_type::const_iterator titer_end = dot.terminals.end();
      for (id_map_type::const_iterator titer = dot.terminals.begin(); titer != titer_end; ++ titer) {
	
	const grammar_node_type& dot_next = grammar_nodes[titer->second];
	
	const bool has_rule = dot_next.edge != hypergraph_type::invalid;
	const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	
	if (has_rule)
	  insert_edge(edge_type(edge.lhs, dot_next, edge, dot_next.edge));
	if (has_next)
	  insert_edge(edge_type(edge.lhs, dot_next, edge));
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
	
	insert_edge(edge_type(lhs, dot_next, edge));
      }
    }

    // comlete passive with actives
    void complete_passive(const edge_type& passive)
    {
      // we will try find actives whose last match with passive's first
      // do we group by passive's lhs?
      
      edge_set_active_type::const_iterator aiter_end = edges_active.end();
      for (edge_set_active_type::const_iterator aiter = edges_active.begin(); aiter != aiter_end; ++ aiter) {
	const edge_type& active = *(*aiter);
	
	const grammar_node_type& dot = *(active.dot);
	id_map_type::const_iterator niter = dot.non_terminals.find(passive.lhs);
	if (niter == dot.non_terminals.end()) continue;
	
	const grammar_node_type& dot_next = grammar_nodes[niter->second];
	
	const bool has_rule = dot_next.edge != hypergraph_type::invalid;
	const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	
	if (has_rule)
	  insert_edge(edge_type(active.lhs, dot_next, active, passive, dot_next.edge));
	if (has_next)
	  insert_edge(edge_type(active.lhs, dot_next, active, passive));
      }
    }
    
    // find passives that can extend the active
    void complete_active(const edge_type& active)
    {
      const grammar_node_type& dot = *(active.dot);
      
      // find passives whose first match with active's last
      
      edge_set_passive_type::const_iterator piter_end = edges_passive.end();
      for (edge_set_passive_type::const_iterator piter = edges_passive.begin(); piter != piter_end; ++ piter) {
	const edge_type& passive = *(*piter);
	
	id_map_type::const_iterator niter = dot.non_terminals.find(passive.lhs);
	if (niter == dot.non_terminals.end()) continue;
	
	const grammar_node_type& dot_next = grammar_nodes[niter->second];
	
	const bool has_rule = dot_next.edge != hypergraph_type::invalid;
	const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	
	if (has_rule)
	  insert_edge(edge_type(active.lhs, dot_next, active, passive, dot_next.edge));
	if (has_next)
	  insert_edge(edge_type(active.lhs, dot_next, active, passive));
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
    
    void connect_edge(const edge_type& edge, const hypergraph_type& source, hypergraph_type& target)
    {
      if (edges_unique.find(&edge) == edges_unique.end()) {
	edges_unique.insert(&edge);
	
	if (edge.is_passive())
	  edges_passive.push_back(&edge);
	else
	  edges_active.push_back(&edge);
	
	agenda_finishing.push_back(&edge);
      }
      
      if (! edge.is_passive()) return;
      
      // Edge is completed or scanned...
      // we will uncover structure...
            
      //
      // we can access rule by source.edges[edge].rule
      //
      
      // we will traverse back until we reach predicted item
      // collect the instance of "passive" item which corresponds to non-terminal!
      std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails;
      
      const edge_type* curr = &edge;
      while (curr && ! curr->is_predicted()) {
	if (curr->passive) {
	  non_terminal_node_set_type::iterator niter = non_terminal_nodes.find(curr->passive);
	  if (niter == non_terminal_nodes.end())
	    throw std::runtime_error("no node?");
	  
	  tails.push_back(niter->second);
	}
	curr = curr->active;
      }
      
      std::reverse(tails.begin(), tails.end());
      
      hypergraph_type::id_type head_id;
      if (edge.lhs == goal_symbol && edge.is_passive()) {
	if (! target.is_valid())
	  target.goal = target.add_node().id;
	
	head_id = target.goal;
      } else {
	non_terminal_node_set_type::iterator niter = non_terminal_nodes.find(&edge);
	if (niter == non_terminal_nodes.end())
	  niter = non_terminal_nodes.insert(std::make_pair(&edge, target.add_node().id)).first;
	
	head_id = niter->second;
      }
      
      hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
      edge_new.rule = source.edges[edge.edge].rule;
      edge_new.features = source.edges[edge.edge].features;
      
      target.connect_edge(edge_new.id, head_id);

      //std::cerr << "connected: " << head_id << " rule: " << *edge_new.rule;
      //std::cerr << " tails: ";
      //std::copy(tails.begin(), tails.end(), std::ostream_iterator<hypergraph_type::id_type>(std::cerr, " "));
      //std::cerr << std::endl;
    }


    void initialize_grammar(const hypergraph_type& source)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;

      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;
      typedef std::vector<edge_set_type, std::allocator<edge_set_type> > node_map_type;
      
      
      edges.clear();
      
      traversals.clear();
      non_terminal_nodes.clear();

      agenda_finishing.clear();
      agenda_exploration.clear();
      
      edges_unique.clear();
      edges_active.clear();
      edges_passive.clear();

#if 0
      spans.clear();
      spans.reserve(source.nodes.size());
      spans.resize(source.nodes.size());

      cicada::span_node(source, spans);
      
      non_terminal_set_type non_terminals(source.nodes.size());
      for (int id = 0; id != source.nodes.size(); ++ id) {
	if (id == source.goal)
	  non_terminals[id] = source.edges[source.nodes[id].edges.front()].rule->lhs.non_terminal();
	else
	  non_terminals[id] = ('['
			       + source.edges[source.nodes[id].edges.front()].rule->lhs.non_terminal_strip()
			       +'_'+boost::lexical_cast<std::string>(spans[id].first)
			       +'_'+boost::lexical_cast<std::string>(spans[id].second)
			       + ']');
	//std::cerr << "non-terminal: " << non_terminals[id] << std::endl;
      }
#endif

#if 1
      node_map_type out_edges(source.nodes.size());
      {
	hypergraph_type::edge_set_type::const_iterator eiter_end = source.edges.end();
	for (hypergraph_type::edge_set_type::const_iterator eiter = source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = *eiter;
	  
	  hypergraph_type::edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	  for (hypergraph_type::edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	    out_edges[*niter].push_back(edge.id);
	}
      }
      
      // assigne pseudo non-terminals
      non_terminal_set_type non_terminals(source.nodes.size());
      for (int id = source.nodes.size() - 1; id >= 0; -- id)
	if (out_edges[id].empty())
	  non_terminals[id] = source.edges[source.nodes[id].edges.front()].rule->lhs.non_terminal();
	else {
	  std::string nodes;
	  
	  bool found = false;
	  for (int pos = 0; pos != source.edges[out_edges[id].front()].tails.size(); ++ pos) {
	    const int antecedent_id = source.edges[out_edges[id].front()].tails[pos];
	    const symbol_type non_terminal = source.edges[source.nodes[antecedent_id].edges.front()].rule->lhs.non_terminal();
	    
	    if (antecedent_id == id)
	      nodes += '@' + non_terminal.non_terminal_strip();
	    else {
	      if (nodes.empty())
		nodes = non_terminal.non_terminal_strip();
	      else
		nodes += '|' + non_terminal.non_terminal_strip();
	    }
	  }
	  
	  non_terminals[id] = ('['
			       + non_terminals[source.edges[out_edges[id].front()].head].non_terminal_strip()
			       + ':' + nodes
			       + ']');
	}
#endif
      
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
	rule_type::symbol_set_type::const_iterator siter_end = edge.rule->source.end();
	for (rule_type::symbol_set_type::const_iterator siter = edge.rule->source.begin(); siter != siter_end; ++ siter) {
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
	
	if (grammar_nodes[niter->second].edge == hypergraph_type::invalid) {
	  grammar_nodes[niter->second].edge = edge.id;
	  grammar_nodes[niter->second].features = edge.features;
	} else 
	  grammar_nodes[niter->second].features += edge.features;
      }
    }

    
  private:  
    symbol_type           goal_symbol;
    grammar_node_set_type grammar_nodes;

    edge_set_type edges;

    traversal_set_type traversals;
    non_terminal_node_set_type non_terminal_nodes;

    agenda_type agenda_finishing;
    agenda_type agenda_exploration;
    
    edge_set_unique_type  edges_unique;
    edge_set_active_type  edges_active;
    edge_set_passive_type edges_passive;

    span_set_type spans;
  };
  
  inline
  void generate_earley(const HyperGraph& source, HyperGraph& target)
  {
    GenerateEarley generater;
      
    generater(source, target);
  }
  
};

#endif
