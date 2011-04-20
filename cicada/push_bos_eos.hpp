// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PUSH_BOS_EOS__HPP__
#define __CICADA__PUSH_BOS_EOS__HPP__ 1

#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <utility>

#include <cicada/hypergraph.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sort.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/mathop.hpp>
#include <utils/bithack.hpp>


namespace cicada
{
  
  struct PushBosEos
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    typedef hypergraph_type::rule_type rule_type;

    typedef Vocab  vocab_type;
    typedef rule_type::symbol_type     symbol_type;
    typedef rule_type::symbol_set_type symbol_set_type;

    typedef std::vector<id_type, std::allocator<id_type> > node_map_type;
    typedef std::vector<id_type, std::allocator<id_type> > queue_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > phrase_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& graph)
    {
      if (! source.is_valid()) {
	graph.clear();
	return;
      }
      
      graph = source;
      
      
      queue_type queue;
      queue_type queue_new;

      phrase_type phrase;
      
      queue.push_back(graph.goal);
      
      {
	node_map_type node_map_bos(graph.nodes.size(), id_type(-1));
	
	while (! queue.empty()) {
	  queue_new.clear();
	  
	  queue_type::const_iterator qiter_end = queue.end();
	  for (queue_type::const_iterator qiter = queue.begin(); qiter != qiter_end; ++ qiter) {
	    const node_type& node = graph.nodes[*qiter];
	    
	    node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	      edge_type& edge = graph.edges[*eiter];
	      
	      const symbol_type&     lhs = edge.rule->lhs;
	      const symbol_set_type& rhs = edge.rule->rhs;
	      
	      if (rhs.front().is_terminal()) {
		if (rhs.front() != vocab_type::BOS) {
		  phrase.clear();
		  phrase.push_back(vocab_type::BOS);
		  phrase.insert(phrase.end(), rhs.begin(), rhs.end());
		  
		  edge.rule = rule_type::create(rule_type(lhs, phrase.begin(), phrase.end()));
		}
	      } else {
		const int __non_terminal_index = rhs.front().non_terminal_index();
		const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, 0, __non_terminal_index - 1);
	      
		const id_type node_id = edge.tails[antecedent_index];
	      
		if (node_map_bos[node_id] == id_type(-1)) {
		  const id_type node_bos_id = graph.add_node().id;
		
		  node_map_bos[node_id] = node_bos_id;
		
		  queue_new.push_back(node_bos_id);
		
		  // we will create edges from node_id
		
		  node_type::edge_set_type::const_iterator aiter_end = graph.nodes[node_id].edges.end();
		  for (node_type::edge_set_type::const_iterator aiter = graph.nodes[node_id].edges.begin(); aiter != aiter_end; ++ aiter) {
		    const edge_type& edge_antecedent = graph.edges[*aiter];
		  
		    edge_type& edge_new = graph.add_edge(edge_antecedent);
		  
		    graph.connect_edge(edge_new.id, node_bos_id);
		  }
		}
	      
		edge.tails[antecedent_index] = node_map_bos[node_id];
	      }
	    }
	  }
	
	  queue.swap(queue_new);
	  queue_new.clear();
	}
      }
      
      queue.clear();
      queue_new.clear();
      
      queue.push_back(graph.goal);
      
      {
	node_map_type node_map_eos(graph.nodes.size(), id_type(-1));
	
	while (! queue.empty()) {
	  queue_new.clear();
	
	  queue_type::const_iterator qiter_end = queue.end();
	  for (queue_type::const_iterator qiter = queue.begin(); qiter != qiter_end; ++ qiter) {
	    const node_type& node = graph.nodes[*qiter];
	  
	    node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	      edge_type& edge = graph.edges[*eiter];

	      const symbol_type&     lhs = edge.rule->lhs;
	      const symbol_set_type& rhs = edge.rule->rhs;
	    
	      if (rhs.back().is_terminal()) {
		if (rhs.back() != vocab_type::EOS) {
		  phrase.clear();
		  phrase.insert(phrase.end(), rhs.begin(), rhs.end());
		  phrase.push_back(vocab_type::EOS);
		  
		  edge.rule = rule_type::create(rule_type(lhs, phrase.begin(), phrase.end()));
		}
	      } else {
		const int __non_terminal_index = rhs.back().non_terminal_index();
		const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, static_cast<int>(edge.tails.size() - 1), __non_terminal_index - 1);
	      
		const id_type node_id = edge.tails[antecedent_index];
	      
		if (node_map_eos[node_id] == id_type(-1)) {
		  const id_type node_eos_id = graph.add_node().id;
		
		  node_map_eos[node_id] = node_eos_id;
		
		  queue_new.push_back(node_eos_id);
		
		  // we will create edges from node_id
		
		  node_type::edge_set_type::const_iterator aiter_end = graph.nodes[node_id].edges.end();
		  for (node_type::edge_set_type::const_iterator aiter = graph.nodes[node_id].edges.begin(); aiter != aiter_end; ++ aiter) {
		    const edge_type& edge_antecedent = graph.edges[*aiter];
		  
		    edge_type& edge_new = graph.add_edge(edge_antecedent);
		  
		    graph.connect_edge(edge_new.id, node_eos_id);
		  }
		}
	      
		edge.tails[antecedent_index] = node_map_eos[node_id];
	      }
	    }
	  }
	
	  queue.swap(queue_new);
	  queue_new.clear();
	}
      }
      
      graph.topologically_sort();
    }
    
  };
  
  
  inline
  void push_bos_eos(const HyperGraph& source, HyperGraph& target)
  {
    PushBosEos()(source, target);
  }
  
  inline
  void push_bos_eos(HyperGraph& graph)
  {
    HyperGraph x;
    push_bos_eos(graph, x);
    graph.swap(x);
  }
};

#endif
