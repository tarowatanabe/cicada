// -*- mode: c++ -*-

#ifndef __CICADA__COMPOSE_TREE__HPP__
#define __CICADA__COMPOSE_TREE__HPP__ 1

#include <vector>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/tree_grammar.hpp>
#include <cicada/tree_transducer.hpp>
#include <cicada/hypergraph.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/hashmurmur.hpp>

#include <google/dense_hash_map>

//
// hypergraph-to-hypergraph transduction by
//
//@InProceedings{zhang-EtAl:2009:EMNLP1,
//    author    = {Zhang, Hui  and  Zhang, Min  and  Li, Haizhou  and  Tan, Chew Lim},
//    title     = {Fast Translation Rule Matching for Syntax-based Statistical Machine Translation},
//    booktitle = {Proceedings of the 2009 Conference on Empirical Methods in Natural Language Processing},
//    month     = {August},
//    year      = {2009},
//    address   = {Singapore},
//    publisher = {Association for Computational Linguistics},
//    pages     = {1037--1045},
//    url       = {http://www.aclweb.org/anthology/D/D09/D09-1108}
//}
//
// The terminologies used in the above paper is somewhat confusing:
//   input-forest, input-hypergraph, encoded hyper-path etc.
//
// The algorithm computes:
//  for each node, try matching with tree fragments
//  when matched, the match is represented by a set of hypergraph-node-id of input-hypergraph
//  and a set of matched rules.
//  We can uncover output hypergraph by enumerating matched rules.
//  Book-keep the matched input-hypergraph in a chart so that we can construct translational packed forest.
// 
// TODO
// support phrasal/synchronous-CFG style rules:
//   We will rune CKY and try match syntactically constrained phrases indicated by source spans...

namespace cicada
{
  
  struct ComposeTree
  {
    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;

    typedef TreeGrammar    grammar_type;
    typedef TreeTransducer transducer_type;
    typedef HyperGraph     hypergraph_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    
    typedef transducer_type::rule_pair_set_type tree_rule_pair_set_type;
    typedef transducer_type::rule_pair_type     tree_rule_pair_type;
    typedef transducer_type::rule_type          tree_rule_type;
    typedef transducer_type::rule_ptr_type      tree_rule_ptr_type;
    
    typedef std::vector<transducer_type::edge_id_type, std::allocator<transducer_type::edge_id_type> > edge_set_type;
    typedef std::vector<transducer_type::id_type, std::allocator<transducer_type::id_type> >           node_set_type;
    
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > frontier_type;
    typedef std::vector<frontier_type, std::allocator<frontier_type> > frontier_set_type;
    
    typedef std::pair<frontier_type, transducer_type::id_type> frontier_pair_type;
    typedef std::deque<frontier_pair_type, std::allocator<frontier_pair_type> > queue_type;

    struct NodeMap
    {
      typedef google::dense_hash_map<symbol_type, hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<symbol_type> > node_map_type;

      typedef node_map_type::value_type value_type;
      
      typedef node_map_type::iterator       iterator;
      typedef node_map_type::const_iterator const_iterator;

      NodeMap() : node_map() { node_map.set_empty_key(symbol_type()); }
      
      std::pair<iterator, bool> insert(const value_type& x) { return node_map.insert(x); }
      
      const_iterator find(const symbol_type& x) const { return node_map.find(x); }
      iterator find(const symbol_type& x) { return node_map.find(x); }
      
      const_iterator begin() const { return node_map.begin(); }
      iterator begin() { return node_map.begin(); }

      const_iterator end() const { return node_map.end(); }
      iterator end() { return node_map.end(); }

      node_map_type node_map;
    };

    typedef NodeMap node_map_type;
    typedef std::vector<node_map_type, std::allocator<node_map_type> > node_map_set_type;
    
    ComposeTree(const symbol_type& __goal, const grammar_type& __grammar, const bool __yield_source)
      : goal(__goal), grammar(__grammar), yield_source(__yield_source) 
    { }
    
    void operator()(const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      graph_out.clear();
      node_map.clear();
      
      if (! graph_in.is_valid()) return;

      node_map.reserve(graph_in.nodes.size());
      node_map.resize(graph_in.nodes.size());
      
      // bottom-up topological order
      for (size_t id = 0; id != graph_in.nodes.size(); ++ id)
	match_tree(id, graph_in, graph_out);
      
      node_map_type::const_iterator niter = node_map[graph_in.goal].find(goal.non_terminal());
      if (niter != node_map[graph_in.goal].end())
	graph_out.goal = niter->second;

      node_map.clear();
    }
    
    void match_tree(const int id, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      if (graph_in.nodes[id].edges.empty()) return;

      queue_type queue;
      
      for (size_t grammar_id = 0; grammar_id != grammar.size(); ++ grammar_id) {
	const transducer_type& transducer = grammar[grammar_id];
	
	const transducer_type::edge_id_type edge_id = transducer.edge(graph_in.edges[graph_in.nodes[id].edges.front()].rule->lhs);
	
	if (edge_id == transducer_type::edge_id_type(-1)) continue;
	
	const transducer_type::edge_id_type edge_epsilon = transducer.edge(vocab_type::EPSILON);
	const transducer_type::edge_id_type edge_none    = transducer.edge(vocab_type::NONE);
	
	transducer_type::id_type node = transducer.next(transducer.root(), edge_id);
	if (node == transducer.root()) continue;
	
	node = transducer.next(node, edge_none);
	if (node == transducer.root()) continue;
	
	queue.clear();
	queue.push_back(std::make_pair(frontier_type(1, id), node));
	
	while (! queue.empty()) {
	  const frontier_pair_type& frontier = queue.front();
	  
	  frontier_set_type frontiers(1, frontier_type());
	  frontier_set_type frontiers_next;
	  
	  node_set_type nodes(1, frontier.second);
	  node_set_type nodes_next;

	  edge_set_type edges;
	  
	  frontier_type::const_iterator niter_end = frontier.first.end();
	  for (frontier_type::const_iterator niter = frontier.first.begin(); niter != niter_end; ++ niter) {
	    
	    frontiers_next.clear();
	    nodes_next.clear();

	    edges.clear();
	    hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = graph_in.nodes[*niter].edges.end();
	    for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = graph_in.nodes[*niter].edges.begin(); eiter != eiter_end; ++ eiter) {
	      const hypergraph_type::edge_type& edge = graph_in.edges[*eiter];
	      
	      edges.push_back(transducer.edge(edge.rule->rhs));
	    }
	    
	    
	    frontier_set_type::const_iterator fiter = frontiers.begin();
	    node_set_type::const_iterator titer_end = nodes.end();
	    for (node_set_type::const_iterator titer = nodes.begin(); titer != titer_end; ++ titer, ++ fiter) {
	      const transducer_type::size_type node_epsilon = transducer.next(edge_epsilon, *titer);
	      if (node_epsilon != transducer.root()) {
		frontier_type frontier(*fiter);
		frontier.push_back(*niter);
		
		frontiers_next.push_back(frontier);
		nodes_next.push_back(node_epsilon);
	      }

	      edge_set_type::const_iterator eiter_begin = edges.begin();
	      edge_set_type::const_iterator eiter_end = edges.end();
	      for (edge_set_type::const_iterator eiter = eiter_begin; eiter != eiter_end; ++ eiter)
		if (*eiter != transducer_type::edge_id_type(-1)) {
		  const transducer_type::edge_id_type& edge_id = *eiter;
		  
		  const transducer_type::size_type node_edge = transducer.next(edge_id, *titer);
		  if (node_edge != transducer.root()) {
		    const hypergraph_type::edge_type& edge = graph_in.edges[graph_in.nodes[*niter].edges[eiter - eiter_begin]];
		    
		    frontier_type frontier(*fiter);
		    frontier.insert(frontier.end(), edge.tails.begin(), edge.tails.end());
		    
		    frontiers_next.push_back(frontier);
		    nodes_next.push_back(node_edge);
		  }
		}
	    }
	    
	    frontiers.swap(frontiers_next);
	    nodes.swap(nodes_next);
	    
	    frontiers_next.clear();
	    nodes_next.clear();
	  }
	  
	  // frontiers and nodes contain new frontier!
	  // in addition, we need to traverse transducer.next() with edge_none!
	  
	  frontier_set_type::const_iterator fiter = frontiers.begin();
	  node_set_type::const_iterator titer_end = nodes.end();
	  for (node_set_type::const_iterator titer = nodes.begin(); titer != titer_end; ++ titer, ++ fiter) {
	    const transducer_type::size_type node_none = transducer.next(edge_none, *titer);
	    if (node_none == transducer.root()) continue;
	    
	    const transducer_type::rule_pair_set_type& rules = transducer.rules(node_none);
	    
	    if (! rules.empty()) {
	      // try match with rules with *fiter == frontier-nodes and generate graph_out!
	      
	      transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	      for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
		apply_rule(*riter, id, *fiter, graph_in, graph_out);
	    }
	    
	    queue.push_back(std::make_pair(*fiter, node_none));
	  }
	  
	  queue.pop_front();
	}
      }
    }
    
    
    void apply_rule(const tree_rule_pair_type& rule_pair,
		    const hypergraph_type::id_type root_in,
		    const frontier_type& frontiers,
		    const hypergraph_type& graph_in,
		    hypergraph_type& graph_out)
    {
      // apply rule for the frontiers with root at graph_in and generate graph_out...
      // frontiers are ordered wrt source side rule's frontier nodes (in pre-traversal order),## the same as post-traversal order...
      //
      // collect frontiers...
      
      //
      // construct graph_out in pre-order...
      //

      const tree_rule_type rule = (yield_source ? *rule_pair.source : *rule_pair.target);
      
      std::pair<node_map_type::iterator, bool> result = node_map[root_in].insert(std::make_pair(rule.label.non_terminal(), 0));
      if (result.second)
	result.first->second = graph_out.add_node().id;
      
      const hypergraph_type::id_type edge_id = construct_graph(rule, result.first->second, frontiers, graph_in, graph_out);
      
      graph_out.edges[edge_id].features += rule_pair.features;
    }
    
    hypergraph_type::id_type construct_graph(const tree_rule_type& rule,
					     const hypergraph_type::id_type root,
					     const frontier_type& frontiers,
					     const hypergraph_type& graph_in,
					     hypergraph_type& graph_out)
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
	    const int non_terminal_index = aiter->label.non_terminal_index();
	    
	    if (non_terminal_index == 0)
	      throw std::runtime_error("invalid non-terminal index");
	    if (non_terminal_index - 1 >= static_cast<int>(frontiers.size()))
	      throw std::runtime_error("non-terminal index exceeds frontier size");
	    
	    const hypergraph_type::id_type node = frontiers[non_terminal_index - 1];

	    std::pair<node_map_type::iterator, bool> result = node_map[node].insert(std::make_pair(aiter->label.non_terminal(), 0));
	    if (result.second)
	      result.first->second = graph_out.add_node().id;
	    
	    tails.push_back(result.first->second);
	  } else
	    tails.push_back(graph_out.add_node().id);
	  
	  rhs.push_back(aiter->label.non_terminal(pos));
	  
	  ++ pos;
	} else
	  rhs.push_back(aiter->label);
      
      const hypergraph_type::id_type edge_id = graph_out.add_edge(tails.begin(), tails.end()).id;
      graph_out.edges[edge_id].rule.reset(new rule_type(rule.label, rule_type::symbol_set_type(rhs.begin(), rhs.end())));
      graph_out.connect_edge(edge_id, root);
      
      tails_type::const_iterator titer = tails.begin();
      for (tree_rule_type::const_iterator aiter = rule.begin(); aiter != aiter_end; ++ aiter) {
	if (aiter->label.is_non_terminal() && ! aiter->antecedents.empty()) {
	  construct_graph(*aiter, *titer, frontiers, graph_in, graph_out);
	  ++ titer;
	}
      }
      
      return edge_id;
    }

    node_map_set_type node_map;
    
    symbol_type goal;
    const grammar_type& grammar;
    const bool yield_source;
  };
  
  
  inline
  void compose_tree(const Symbol& goal, const TreeGrammar& grammar, const HyperGraph& graph_in, HyperGraph& graph_out, const bool yield_source=false)
  {
    ComposeTree __composer(goal, grammar, yield_source);
    __composer(graph_in, graph_out);
  }
};

#endif
