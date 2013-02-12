// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__VERIFY__HPP__
#define __CICADA__VERIFY__HPP__ 1

#include <vector>
#include <stdexcept>

#include <cicada/hypergraph.hpp>

#include <utils/bithack.hpp>

namespace cicada
{
  struct Verify
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    typedef Vocab      vocab_type;
    
    typedef hypergraph_type::rule_type    rule_type;
    
    typedef rule_type::symbol_type     symbol_type;
    typedef rule_type::symbol_set_type symbol_set_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > label_set_type;
    typedef std::vector<bool, std::allocator<bool> > visited_set_type;
    
    void operator()(const hypergraph_type& graph)
    {
      // verify goal...
      if (! graph.is_valid())
	throw std::runtime_error("invalid hypergraph: no goal?");

      label_set_type   labels(graph.nodes.size());
      visited_set_type visited(graph.edges.size(), false);

      // verify nodes...
      size_type node_pos = 0;
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter, ++ node_pos) {
	const hypergraph_type::node_type& node = *niter;
	
	if (node_pos != node.id)
	  throw std::runtime_error("node id does not match?");
	
	if (node.edges.empty())
	  throw std::runtime_error("no hyperedges in a node?");
	
	// verify antecedents...
	bool initial = true;
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  if (*eiter >= graph.edges.size())
	    throw std::runtime_error("invalid hyperedge id");

	  if (visited[*eiter])
	    throw std::runtime_error("multiple instances of hyperedge?");
	  
	  visited[*eiter] = true;
	  
	  if (graph.edges[*eiter].head != node_pos)
	    throw std::runtime_error("invalid head in a hyperedge?");
	  
	  if (! graph.edges[*eiter].rule)
	    throw std::runtime_error("no rules associated with a hyperedge?");
	  
	  if (initial)
	    labels[node_pos] = graph.edges[*eiter].rule->lhs;
	  else if (labels[node_pos] != graph.edges[*eiter].rule->lhs)
	    throw std::runtime_error("invalid node-label assigned");
	  
	  initial = false;
	}
      }

      // verify edge
      size_type edge_pos = 0;
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter, ++ edge_pos) {
	const hypergraph_type::edge_type& edge = *eiter;
	
	if (edge_pos != edge.id)
	  throw std::runtime_error("edge id does not match?");
	
	if (! edge.rule)
	  throw std::runtime_error("no rules associated with a hyperedge?");

	if (! visited[edge_pos])
	  throw std::runtime_error("isolated hyperedge?");
	
	// verify tails...
	hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	  if (*titer >= graph.nodes.size())
	    throw std::runtime_error("invalid tail id?");
	
	int non_terminal_pos = 0;
	symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	for (symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
	  if (riter->is_non_terminal()) {
	    const int __non_terminal_index = riter->non_terminal_index();
	    const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	    ++ non_terminal_pos;
	    
	    if (antecedent_index >= static_cast<int>(edge.tails.size()))
	      throw std::runtime_error("invalid tails and rule's rhs?");
	    
	    if (riter->non_terminal() != labels[edge.tails[antecedent_index]])
	      throw std::runtime_error("invalid label in rhs?");
	  }
	
	if (non_terminal_pos != static_cast<int>(edge.tails.size()))
	  throw std::runtime_error("rhs size and tails size does not match?");
      }
    }
  };

  inline
  void verify(HyperGraph& graph)
  {
    Verify __verify;
    __verify(graph);
  }

};

#endif
