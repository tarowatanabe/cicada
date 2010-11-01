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
    
    ComposeTree(const grammar_type& __grammar)
      : grammar(__grammar) {}
    
    void operator()(const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      graph_out.clear();

      if (! graph_in.is_valid()) return;
      
      // bottom-up topological order
      for (size_t id = 0; id != graph_in.size(); ++ id)
	match_tree(id, graph_in, graph_out);
      
      node_map_type::const_iterator niter = node_map.find(graph_in.goal);
      if (niter != node_map.end())
	graph_out.goal = niter->second;
    }
    
    void match_tree(const int id, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      if (graph_in.nodes[id].edges.empty()) return;
      
      for (size_t grammar_id = 0; grammar_id != grammar.size(); ++ id) {
	const transducer_type& transducer = grammar[grammar_id];
	
	const transducer_type::edge_id_type edge_id = transducer.edge(graph_in.edges[graph_in.nodes[id].edges.front()].rule->lhs);
	
	if (edge_id == transducer_type::edge_id_type(-1)) continue;
	
	const transducer_type::edge_id_type edge_epsilon = transducer.edge(vocab_type::EPSILON);
	const transducer_type::edge_id_type edge_none    = transducer.edge(vocab_type::NONE);
	
	transducer_type::id_type node = transducer.next(transducer.root(), edge_id);
	if (node == transducer.root()) continue;
	
	node = transducer.next(node, edge_none);
	if (node == transducer.root()) continue;
	
	queue.push_back(frontier_type(node_set_type(id), node));
	
	while (! queue.empty()) {
	  const frontier_type& frontier = queue.front();
	  
	  frontier_set_type frontiers(1, node_set_type());
	  frontier_set_type frontiers_next;
	  
	  transducer_node_set_type transducer_nodes(1, frontier.second);
	  transducer_node_set_tyoe transducer_nodes_next;

	  edge_set_type edges;
	  
	  node_set_type::const_iterator niter_end = frontier.first.end();
	  for (node_set_type::const_iterator niter = frontier.first.begin(); niter != niter_end; ++ niter) {
	    
	    frontiers_next.clear();
	    transducer_nodes_next.clear();

	    edges.clear();
	    hypergraph_type::node_type::edge_set_type eiter_end = graph_in.node[*niter].edges.end();
	    for (hypergraph_type::node_type::edge_set_type eiter = graph_in.node[*niter].edges.begin(); eiter != eiter_end; ++ eiter) {
	      const hypergraph_type::edge_tyep& edge = graph_in.edges[*eiter];
	      
	      edges.push_back(transducer.edge(edge.rule->rhs));
	    }
	    
	    
	    node_set_type::const_iterator fiter = frontiers.begin();
	    transducer_node_set_type::const_iterator titer_end = transducer_nodes.end();
	    for (transducer_node_set_type::const_iterator titer = transducer_nodes.begin(); titer != titer_end; ++ titer, ++ fiter) {
	      const transducer_type::size_type node_epsilon = transducer.next(edge_epsilon, *titer);
	      if (node_epsilon != transducer.root()) {
		frontier_type frontier(*fiter);
		frontier.push_back(*niter);
		
		frontiers_next.push_back(frontier);
		transducer_nodes_next.push_back(node_epsilon);
	      }

	      edge_set_type::const_iterator eiter_begin = edges.begin();
	      edge_set_type::const_iterator eiter_end = edges.end();
	      for (edge_set_type::const_iterator eiter = eiter_begin; eiter != eiter_end; ++ eiter)
		if (*eiter != transducer_type::edge_id_type(-1)) {
		  const transducer_type::edge_id_type& edge_id = *eiter;
		  
		  const transducer_type::size_type node_edge = transducer.next(edge_id, *titer);
		  if (node_edge != transducer.root()) {
		    const hypergraph_type::edge_type& edge = graph_in.edges[graph_in.node[*niter].edges[eiter - eiter_begin]];
		    
		    frontier_type frontier(*fiter);
		    frontier.insert(frontier.end(), edge.tails.begin(), edge.tails.end());
		    frontiers_next.push_back(frontier);
		    transducer_nodes_next.push_back(node_edge);
		  }
		}
	    }
	    
	    frontiers.swap(frontiers_next);
	    transducer_nodes.swap(transducer_nodes_next);
	    
	    frontiers_next.clear();
	    transducer_nodes_next.clear();
	  }
	  
	  // frontiers and transducer_nodes contain new frontier!
	  // in addition, we need to traverse transducer.next() with edge_none!
	  
	  node_set_type::const_iterator fiter = frontiers.begin();
	  transducer_node_set_type::const_iterator titer_end = transducer_nodes.end();
	  for (transducer_node_set_type::const_iterator titer = transducer_nodes.begin(); titer != titer_end; ++ titer, ++ fiter) {
	    const transducer_type::size_type node_none = transducer.next(edge_none, *titer);
	    if (node_none == transducer.root()) continue;
	    
	    const transducer_type::rule_pair_set_tyep& rules = transducer.rules(node_none);
	    
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
    
    
    void apply_rule(const tree_rule_pair_type& rule_pair, const int root, const node_set_type& frontiers, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      // apply rule for the frontiers with root at graph_in and generate graph_out...
      // frontiers are ordered wrt source side rule's frontier nodes (in pre-traversal order),## the same as post-traversal order...
      //
      // collect frontiers...
      
      //
      // construct graph_out in pre-order...
      //
      
      node_map_type::iterator niter = node_map.find(root);
      if (niter == node_map.end())
	niter->second = node_map.insert(std::make_pair(root, graph_out.add_node().id())).first;
      
      const hypergraph_type::id_type edge_id = construct_graph(*rule_pair.target, niter->second, frontiers, graph_in, graph_out);
      
      graph_out.edges[edge_id].features += rule_pair.features;
    }
    
    hypergraph_type::id_type construct_graph(const tree_rule_type& rule,
					     const hypergraph_type::id_type root,
					     const node_set_type& frontiers,
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
	    
	    const hypergraph_type::id_type node = frontiers[non_terminal_index - 1];
	    
	    node_map_type::iterator niter = node_map.find(node);
	    if (niter == node_map.end())
	      niter->second = node_map.insert(std::make_pair(node, graph_out.add_node().id())).first;
	    
	    tails.push_back(niter->second);
	  } else
	    tails.push_back(graph_out.add_node().id());
	  
	  rhs.push_back(aiter->label.non_terminal(pos));
	  
	  ++ pos;
	} else
	  rhs.push_back(aiter->label);
      
      const hypegraph_type::id_type edge_id = graph_out.add_edge(tails.begin(), tails.end()).id;
      graph_out.edges[edge_id].rule.reset(new rule_type(rule.label, rule_type::symbol_set_type(rhs.begin(), rhs.end())));
      graph_out.connect_edge(edge_id, root);
      
      tails_type::const_iterator titer = tails.begin();
      tree_rule_type::const_iterator aiter_end = rule.end();
      for (tree_rule_type::const_iterator aiter = rule.begin(); aiter != aiter_end; ++ aiter) {
	if (aiter->label.is_non_terminal() && ! aiter->antecedents.empty()) {
	  construct_graph(*aiter, *titer, frontiers, graph_in, graph_out);
	  ++ titer;
	}
      }
      
      return edge_id;
    }
    
    const grammar_type& grammar;
  };
  
  
  inline
  void compose_tree(const Grammar& grammar, const HyperGraph& graph_in, HyperGraph& graph_out)
  {
    ComposeTree __composer(grammar);
    __composer(graph_in, graph_out);
  }
};

#endif
