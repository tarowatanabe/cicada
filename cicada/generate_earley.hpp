// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

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
#include <cicada/span_vector.hpp>
#include <cicada/sentence.hpp>

#include <cicada/inside_outside.hpp>
#include <cicada/semiring.hpp>

#include <google/dense_hash_set>
#include <google/dense_hash_map>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/compact_trie.hpp>
#include <utils/indexed_trie.hpp>
#include <utils/bithack.hpp>

namespace cicada
{
  
  struct GenerateEarley
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;
    
    typedef Sentence   sentence_type;
    typedef Lattice    lattice_type;
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;
    
    GenerateEarley(const int __depth, const int __width)
      : depth(__depth), width(__width)
    {
      items_unique.set_empty_key(0);

      traversals.set_empty_key(traversal_type());

      non_terminal_nodes.set_empty_key(0);
    }
    
    //
    // generate source hypergraph with FST grammar in grammar...!
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

    typedef std::vector<int, std::allocator<int> > depth_set_type;
    
    struct Item
    {
      typedef Item item_type;
      
      symbol_type lhs; // lhs of rule
      const grammar_node_type* dot;  // earley grammar position
      
      // tree depth
      int          depth;
      
      // backpointers
      const item_type* active;
      const item_type* passive;
      
      // backptr to source hypergraph's edge
      hypergraph_type::id_type edge;
      
      // created by predict...
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const int& __depth)
	: lhs(__lhs), dot(&__dot),
	  depth(__depth), 
	  active(0), passive(0),
	  edge(hypergraph_type::invalid) {}
#if 0      
      // created by predict...
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const int& __depth,
	   const item_type& __active)
	: lhs(__lhs), dot(&__dot),
	  depth(__depth), 
	  active(&__active), passive(0),
	  edge(hypergraph_type::invalid) {}
#endif
      // created by scan... we will always have terminal
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const int& __depth, 
	   const item_type& __active,
	   const hypergraph_type::id_type& __edge)
	: lhs(__lhs), dot(&__dot),
	  depth(__depth), 
	  active(&__active), passive(0),
	  edge(__edge) {}
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const int& __depth, 
	   const item_type& __active)
	: lhs(__lhs), dot(&__dot),
	  depth(__depth), 
	  active(&__active), passive(0),
	  edge(hypergraph_type::invalid) {}
      
      // construct by complete
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const int& __depth, 
	   const item_type& __active, const item_type& __passive,
	   const hypergraph_type::id_type& __edge)
	: lhs(__lhs), dot(&__dot),
	  depth(__depth), 
	  active(&__active), passive(&__passive),
	  edge(__edge) {}
      Item(const symbol_type& __lhs, const grammar_node_type& __dot,
	   const int& __depth, 
	   const item_type& __active, const item_type& __passive)
	: lhs(__lhs), dot(&__dot),
	  depth(__depth), 
	  active(&__active), passive(&__passive),
	  edge(hypergraph_type::invalid) {}
      
      // for query...
      Item(const int& __depth) : depth(__depth) {}
      
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
    
    struct item_unique_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const item_type* x) const
      {
	// passive item do not care about what dot we have consumed...
	if (! x)
	  return 0;
	else if (x->is_active())
	  return hasher_type::operator()(x->dot, hasher_type::operator()(x->depth, x->lhs.id()));
	else
	  return hasher_type::operator()(x->depth, x->lhs.id());
      }
    };
    
    struct item_unique_equal_type
    {
      bool operator()(const item_type* x, const item_type* y) const
      {
	// passive item do not care dot...
	return ((x == y) 
		|| (x && y
		    && x->is_active() == y->is_active()
		    && x->lhs == y->lhs
		    && x->depth == y->depth
		    && (x->is_passive() || x->dot == y->dot)));
      }
    };
    
    struct item_active_hash_type : public utils::hashmurmur<size_t>
    {
      size_t operator()(const item_type* x) const
      {
	return utils::hashmurmur<size_t>::operator()(x->depth);
      }
    };
    
    struct item_active_equal_type
    {
      bool operator()(const item_type* x, const item_type* y) const
      {
	return x->depth == y->depth;
      }
    };

    struct item_passive_hash_type : public utils::hashmurmur<size_t>
    {
      size_t operator()(const item_type* x) const
      {
	return utils::hashmurmur<size_t>::operator()(x->depth);
      }
    };
    
    struct item_passive_equal_type
    {
      bool operator()(const item_type* x, const item_type* y) const
      {
	return x->depth == y->depth;
      }
    };

    struct item_node_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const item_type* x) const
      {
	return (x ? hasher_type::operator()(x->depth, x->lhs.id()) : size_t(0));
      }
    };

    struct item_node_equal_type
    {
      bool operator()(const item_type* x, const item_type* y) const
      {
	return ((x == y) || (x && y
			     && x->lhs == y->lhs
			     && x->depth == y->depth));
      }
    };
    
    typedef google::dense_hash_set<const item_type*, item_unique_hash_type, item_unique_equal_type > item_set_unique_type;

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
    typedef google::dense_hash_map<const item_type*, hypergraph_type::id_type, item_node_hash_type, item_node_equal_type > non_terminal_node_set_type;
    
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
      insert_item(item_type(goal_symbol, dot_goal, 0));
      
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
      
      // sort!
      if (target.is_valid())
	target.topologically_sort();
    }

  private:
    
    void scan(const item_type& item)
    {
      const grammar_node_type& dot = *item.dot;

      //
      // scanning implies moving dot on terminals w/o intersecting with transducer...
      //
      
      id_map_type::const_iterator titer_end = dot.terminals.end();
      for (id_map_type::const_iterator titer = dot.terminals.begin(); titer != titer_end; ++ titer) {
	
	const grammar_node_type& dot_next = grammar_nodes[titer->second];
	
	const bool has_rule = dot_next.edge != hypergraph_type::invalid;
	const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	

	// actually, we need to increment item.span.second! but temporary disabled
	if (has_rule)
	  insert_item(item_type(item.lhs, dot_next, item.depth, item, dot_next.edge));
	if (has_next)
	  insert_item(item_type(item.lhs, dot_next, item.depth, item));
      }
    }
    
    void predict(const item_type& item)
    {
      if (item.depth >= max_tree_depth) return;
      
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
	
	insert_item(item_type(lhs, dot_next, item.depth + 1, item));
      }
    }

    // comlete passive with actives
    void complete_passive(const item_type& passive)
    {
      // we will try find actives whose last match with passive's first
      // do we group by passive's lhs?

      const item_type query(passive.depth - 1);
      
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
	  insert_item(item_type(active.lhs, dot_next, active.depth, active, passive, dot_next.edge));
	if (has_next)
	  insert_item(item_type(active.lhs, dot_next, active.depth, active, passive));
      }
    }
    
    // find passives that can extend the active
    void complete_active(const item_type& active)
    {
      const grammar_node_type& dot = *(active.dot);
      
      // find passives whose first match with active's last
      const item_type query(active.depth + 1);
      
      std::pair<item_set_passive_type::const_iterator, item_set_passive_type::const_iterator> result = items_passive.equal_range(&query);
      for (item_set_passive_type::const_iterator piter = result.first; piter != result.second; ++ piter) {
	const item_type& passive = *(*piter);
	
	id_map_type::const_iterator niter = dot.non_terminals.find(passive.lhs);
	if (niter == dot.non_terminals.end()) continue;
	
	const grammar_node_type& dot_next = grammar_nodes[niter->second];
	
	const bool has_rule = dot_next.edge != hypergraph_type::invalid;
	const bool has_next = ! dot_next.terminals.empty() || ! dot_next.non_terminals.empty();
	
	if (has_rule)
	  insert_item(item_type(active.lhs, dot_next, active.depth, active, passive, dot_next.edge));
	if (has_next)
	  insert_item(item_type(active.lhs, dot_next, active.depth, active, passive));
      }
    }
    
    void insert_item(const item_type& item)
    {
      if (item.passive && item.active) {
	if (traversals.find(traversal_type(item.active, item.passive, item.is_active())) != traversals.end())
	  return;
	else
	  traversals.insert(traversal_type(item.active, item.passive, item.is_active()));
      }
      
      items.push_back(item);
      
      agenda_exploration.push_back(&items.back());
    }
    
    void connect_item(const item_type& item, const hypergraph_type& source, hypergraph_type& target)
    {
      if (items_unique.find(&item) == items_unique.end()) {
	items_unique.insert(&item);
	
	if (item.is_passive())
	  items_passive.insert(&item);
	else
	  items_active.insert(&item);
	
	agenda_finishing.push_back(&item);
      }
      
      if (! item.is_passive()) return;
      
      // Item is completed or scanned...
      // we will uncover structure...
            
      //
      // we can access rule by source.edges[edge].rule
      //
      
      // we will traverse back until we reach predicted item
      // collect the instance of "passive" item which corresponds to non-terminal!
      std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails;
      
      const item_type* curr = &item;
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
      if (item.lhs == goal_symbol) {
	if (! target.is_valid())
	  target.goal = target.add_node().id;
	
	head_id = target.goal;
      } else {
	std::pair<non_terminal_node_set_type::iterator, bool> result = non_terminal_nodes.insert(std::make_pair(&item, 0));
	if (result.second) 
	  result.first->second = target.add_node().id;
	
	head_id = result.first->second;
      }
      
      hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
      edge_new.rule = source.edges[item.edge].rule;
      edge_new.features = item.dot->features;
      
      target.connect_edge(edge_new.id, head_id);

      //std::cerr << "connected: " << head_id << " rule: " << *edge_new.rule;
      //std::cerr << " tails: ";
      //std::copy(tails.begin(), tails.end(), std::ostream_iterator<hypergraph_type::id_type>(std::cerr, " "));
      //std::cerr << std::endl;
    }

    struct length_function
    {
      typedef cicada::Vocab vocab_type;
      typedef cicada::Rule rule_type;
      typedef cicada::semiring::Tropical<int> value_type;
      
      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	int length = 0;
	rule_type::symbol_set_type::const_iterator siter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator siter = edge.rule->rhs.begin(); siter != siter_end; ++ siter)
	  length += (*siter != vocab_type::EPSILON && siter->is_terminal());
	
	// since we will "max" at operator+, we will collect negative length
	return cicada::semiring::traits<value_type>::exp(length);
      }
    };


    void initialize_grammar(const hypergraph_type& source)
    {
      typedef std::vector<std::string, std::allocator<std::string> > node_label_set_type;
      typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
      
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_set_type;
      
      typedef std::vector<edge_set_type, std::allocator<edge_set_type> > node_map_type;
      
      items.clear();
      
      traversals.clear();
      non_terminal_nodes.clear();

      agenda_finishing.clear();
      agenda_exploration.clear();
      
      items_unique.clear();
      items_active.clear();
      items_passive.clear();

      int sentence_length = 0;
      {
	std::vector<length_function::value_type, std::allocator<length_function::value_type> > lengths(source.nodes.size());
	
	cicada::inside(source, lengths, length_function());
	
	sentence_length = log(lengths.back());
      }
      
      // max sentence length is 1.5 of max-sentence-length in the hyperraph
      max_sentence_length = sentence_length + (sentence_length >> 1);

      // compute out-edges...
      node_map_type out_edges(source.nodes.size());
      {
	hypergraph_type::edge_set_type::const_iterator eiter_end = source.edges.end();
	for (hypergraph_type::edge_set_type::const_iterator eiter = source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = *eiter;

	  hypergraph_type::edge_type::node_set_type::const_iterator titer_begin = edge.tails.begin();
	  hypergraph_type::edge_type::node_set_type::const_iterator titer_end   = edge.tails.end();
	  
	  for (hypergraph_type::edge_type::node_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	    out_edges[*titer].push_back(edge.id);
	}
      }
      
      // compute node labels and depth
      // how to compute left/right items?
      non_terminal_set_type non_terminals(source.nodes.size());
      node_label_set_type labels(source.nodes.size());
      depth_set_type depths(source.nodes.size(), 0);

      max_tree_depth = 0;
      
      for (int id = source.nodes.size() - 1; id >= 0; -- id) {
	std::string& label = labels[id];
	
	if (out_edges[id].empty())
	  label = source.edges[source.nodes[id].edges.front()].rule->lhs.non_terminal_strip();
	else {
	  const hypergraph_type::edge_type& edge_parent = source.edges[out_edges[id].front()];
	  const hypergraph_type::id_type parent_id = edge_parent.head;
	  
	  depths[id] = depths[parent_id] + 1;
	  max_tree_depth = utils::bithack::max(max_tree_depth, depths[id]);
	  
	  if (width == 0)
	    label = source.edges[source.nodes[id].edges.front()].rule->lhs.non_terminal().non_terminal_strip();
	  else if (width < 0) {
	    hypergraph_type::edge_type::node_set_type::const_iterator titer_begin = edge_parent.tails.begin();
	    hypergraph_type::edge_type::node_set_type::const_iterator titer_end   = edge_parent.tails.end();
	    for (hypergraph_type::edge_type::node_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	      const int antecedent_id = *titer;
	      const symbol_type non_terminal = source.edges[source.nodes[antecedent_id].edges.front()].rule->lhs.non_terminal();
	      
	      if (antecedent_id == id)
		label += '@' + non_terminal.non_terminal_strip();
	      else if (label.empty())
		label = non_terminal.non_terminal_strip();
	      else
		label += '|' + non_terminal.non_terminal_strip();
	    }
	  } else {
	    hypergraph_type::edge_type::node_set_type::const_iterator titer_begin = edge_parent.tails.begin();
	    hypergraph_type::edge_type::node_set_type::const_iterator titer_end   = edge_parent.tails.end();
	    hypergraph_type::edge_type::node_set_type::const_iterator titer_curr = titer_end;
	    
	    for (hypergraph_type::edge_type::node_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	      if (static_cast<int>(*titer) == id) {
		titer_curr = titer;
		break;
	      }
	    
	    if (titer_curr == titer_end)
	      throw std::runtime_error("no antecedent?");
	    
	    hypergraph_type::edge_type::node_set_type::const_iterator titer_first = std::max(titer_begin, titer_curr - width);
	    for (hypergraph_type::edge_type::node_set_type::const_iterator titer = titer_first; titer != titer_curr; ++ titer) {
	      const symbol_type non_terminal = source.edges[source.nodes[*titer].edges.front()].rule->lhs.non_terminal();
	      
	      if (label.empty())
		label = non_terminal.non_terminal_strip();
	      else
		label += '|' + non_terminal.non_terminal_strip();
	    }
	    
	    label += '@' + source.edges[source.nodes[*titer_curr].edges.front()].rule->lhs.non_terminal().non_terminal_strip();
	    
	    hypergraph_type::edge_type::node_set_type::const_iterator titer_last = std::min(titer_curr + 1 + width, titer_end);
	    for (hypergraph_type::edge_type::node_set_type::const_iterator titer = titer_curr + 1; titer != titer_last; ++ titer) {
	      const symbol_type non_terminal = source.edges[source.nodes[*titer].edges.front()].rule->lhs.non_terminal();
	      
	      label += '|' + non_terminal.non_terminal_strip();
	    }
	  }
	  
	  if (depth != 1) {
	    label = labels[source.edges[out_edges[id].front()].head] + ';' + label;
	    // always strip-off the first non-terminal up-until ';'
	    if (depth > 0 && depths[id] > depth) {
	      std::string::size_type pos = label.find(';', 1);
	      if (pos != std::string::npos)
		label = label.substr(pos + 1);
	    }
	  }
	}
	
	non_terminals[id] = '[' + labels[id] + ']';
	
	//std::cerr << "non-terminal: " << non_terminals[id] << std::endl;
      }
      
      // max-tree-depth is the 1.5 of the original tree
      max_tree_depth = max_tree_depth + (max_tree_depth >> 1);
      
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
	    const int __non_terminal_index = siter->non_terminal_index();
	    const int pos = utils::bithack::branch(__non_terminal_index <= 0, tail_pos, __non_terminal_index - 1);
	    ++ tail_pos;
	    
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
	
	if (grammar_nodes[niter->second].edge == hypergraph_type::invalid) {
	  grammar_nodes[niter->second].edge = edge.id;
	  grammar_nodes[niter->second].features = edge.features;
	} else
	  grammar_nodes[niter->second].features += edge.features;
      }
    }

    
  private:  
    int depth;
    int width;

    int max_sentence_length;
    int max_tree_depth;

    symbol_type           goal_symbol;
    grammar_node_set_type grammar_nodes;

    item_set_type items;

    traversal_set_type traversals;
    non_terminal_node_set_type non_terminal_nodes;

    agenda_type agenda_finishing;
    agenda_type agenda_exploration;
    
    item_set_unique_type  items_unique;
    item_set_active_type  items_active;
    item_set_passive_type items_passive;
  };
  
  inline
  void generate_earley(const HyperGraph& source, HyperGraph& target, const int depth=0, const int width=0)
  {
    GenerateEarley generater(depth, width);
      
    generater(source, target);
  }
  
};

#endif
