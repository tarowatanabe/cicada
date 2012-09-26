// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PUSH_HEAD__HPP__
#define __CICADA__PUSH_HEAD__HPP__ 1

//
// push heads... actually, we will assing "head" attribute...
//

#include <string>
#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <utils/bithack.hpp>

namespace cicada
{
  struct PushHead
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    typedef Vocab      vocab_type;
    typedef Symbol     symbol_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;
    
    typedef hypergraph_type::feature_set_type::feature_type     feature_type;
    typedef hypergraph_type::attribute_set_type::attribute_type attribute_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > head_set_type;
    
    PushHead() : attr_head("head") {}
    
    const attribute_type attr_head;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;
      
      //
      // starting from root, then propagate toward the leafs...
      // if we encountered multiple terminals, use the left-most one...
      //

      head_set_type heads(target.nodes.size(), vocab_type::EPSILON);
      
      hypergraph_type::node_set_type::const_reverse_iterator niter_end = target.nodes.rend();
      for (hypergraph_type::node_set_type::const_reverse_iterator niter = target.nodes.rbegin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  hypergraph_type::edge_type& edge = target.edges[*eiter];
	  
	  // assign head!
	  edge.attributes[attr_head] = static_cast<const std::string&>(heads[node.id]);
	  
	  symbol_type head = heads[node.id];
	  
	  rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
	    if (riter->is_terminal()) {
	      head = *riter;
	      break;
	    }
	  
	  if (head == vocab_type::EPSILON) continue;
	  
	  edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	  for (edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer) {
	    if (heads[*titer] == vocab_type::EPSILON)
	      heads[*titer] = head;
	    else if (heads[*titer] != head)
	      throw std::runtime_error("invalid head pushing!");
	  }
	}
      }
    }    
  };
  
  inline
  void push_head(const HyperGraph& source, HyperGraph& target)
  {
    PushHead()(source, target);
  }
  
  inline
  void push_head(HyperGraph& graph)
  {
    HyperGraph x;
    push_head(graph, x);
    graph.swap(x);
  }

};

#endif
