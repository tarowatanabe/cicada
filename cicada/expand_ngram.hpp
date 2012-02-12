// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EXPAND_NGRAM__HPP__
#define __CICADA__EXPAND_NGRAM__HPP__ 1

#include <vector>
#include <utility>

#include <cicada/hypergraph.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/small_vector.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/bithack.hpp>

namespace cicada
{
   
  class ExpandNGram
  {
  public:
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type   id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef rule_type::vocab_type      vocab_type;
    typedef rule_type::symbol_type     symbol_type;
    typedef rule_type::symbol_set_type symbol_set_type;
    
    typedef symbol_type word_type;
    
    typedef std::vector<word_type, std::allocator<word_type> > buffer_type;
    
    typedef symbol_set_type context_type;
    
    typedef std::pair<context_type, context_type> context_pair_type;
    typedef std::vector<context_pair_type, std::allocator<context_pair_type> > context_pair_set_type;

    typedef std::vector<id_type, std::allocator<id_type> > node_set_type;
    typedef std::vector<node_set_type, std::allocator<node_set_type> > node_map_type;

    typedef context_pair_type state_type;
    
    struct state_hash_type : public utils::hashmurmur<size_t>
    {
      size_t operator()(const state_type& x) const
      {
	typedef utils::hashmurmur<size_t> hasher_type;
	
	return hasher_type::operator()(x.first.begin(), x.first.end(), hasher_type::operator()(x.second.begin(), x.second.end(), 0));
      }
    };

#ifdef HAVE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_map<state_type, id_type, state_hash_type, std::equal_to<state_type>, std::allocator<std::pair<const state_type, id_type> > > state_set_type;
#else
    typedef sgi::hash_map<state_type, id_type, state_hash_type, std::equal_to<state_type>, std::allocator<std::pair<const state_type, id_type > > > state_set_type;
#endif
    
    typedef utils::small_vector<int, std::allocator<int> > index_set_type;

    
  public:
    ExpandNGram(const int __order)
      : order(__order) {}
    
    void operator()(const hypergraph_type& graph, hypergraph_type& graph_out)
    {
      graph_out.clear();
      
      node_map.clear();
      node_map.reserve(graph.nodes.size());
      node_map.resize(graph.nodes.size());
      
      contexts.clear();
      contexts.reserve(graph.nodes.size() * 128);
      
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
	const node_type& node = *niter;
	const bool is_goal = (node.id == graph.goal);
	
	state_buf.clear();
	
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const edge_type& edge = graph.edges[*eiter];
	  
	  index_set_type j_ends(edge.tails.size(), 0);
	  index_set_type j(edge.tails.size(), 0);
	  
	  for (size_t i = 0; i != edge.tails.size(); ++ i)
	    j_ends[i] = node_map[edge.tails[i]].size();
	  
	  edge_type::node_set_type tails(edge.tails.size());
	  
	  for (;;) {
	    // current tails...
	    for (size_t i = 0; i != edge.tails.size(); ++ i)
	      tails[i] = node_map[edge.tails[i]][j[i]];

	    edge_type& edge_new = graph_out.add_edge(tails.begin(), tails.end());
	    edge_new.head = edge.head;
	    edge_new.rule = edge.rule;
	    edge_new.features   = edge.features;
	    edge_new.attributes = edge.attributes;
	    
	    // apply various ngram cconetxt...
	    const state_type state = apply(edge, tails);
	    
	    if (is_goal) {
	      if (graph_out.goal == hypergraph_type::invalid) {
		graph_out.goal = graph_out.add_node().id;
		contexts.push_back(state);
	      } 
	      
	      graph_out.connect_edge(edge_new.id, graph_out.nodes[graph_out.goal].id);
	    } else {
	      typedef std::pair<state_set_type::iterator, bool > result_type;
	      
	      result_type result = state_buf.insert(std::make_pair(state, 0));
	      if (result.second) {
		result.first->second = graph_out.add_node().id;
		
		contexts.push_back(state);
		node_map[edge.head].push_back(result.first->second);
	      }
	      
	      graph_out.connect_edge(edge_new.id, result.first->second);
	    }
	    
	    // proceed to the next j
	    size_t index = 0;
	    for (/**/; index != edge.tails.size(); ++ index) {
	      ++ j[index];
	      if (j[index] < j_ends[index]) break;
	      j[index] = 0;
	    }
	    // finished!
	    if (index == edge.tails.size()) break;
	  }
	}
      }
      
      // topologically sort...
      if (graph_out.is_valid())
	graph_out.topologically_sort();
    }

    template <typename Tails>
    state_type apply(const edge_type& edge, const Tails& tails)
    {
      const context_type& context = edge.rule->rhs;

      const int context_size = order - 1;
      
      buffer.clear();
      
      if (tails.empty()) {
	context_type::const_iterator citer_end = context.end();
	for (context_type::const_iterator citer = context.begin(); citer != citer_end; ++ citer)
	  if (*citer != vocab_type::EPSILON)
	    buffer.push_back(*citer);
	
	const state_type state(static_cast<int>(buffer.size()) <= context_size
			       ? std::make_pair(context_type(buffer.begin(), buffer.end()),
						context_type())
			       : std::make_pair(context_type(buffer.begin(), buffer.begin() + context_size),
						context_type(buffer.end() - context_size, buffer.end())));
	
	return state;
      } else {
	buffer.reserve(context.size() + tails.size() * order * 2);
	
	int star_first = -1;
	int star_last  = -1;
	
	int non_terminal_pos = 0;
	context_type::const_iterator citer_end = context.end();
	for (context_type::const_iterator citer = context.begin(); citer != citer_end; ++ citer) {
	  if (citer->is_non_terminal()) {
	    const int __non_terminal_index = citer->non_terminal_index();
	    const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	    ++ non_terminal_pos;
	    
	    const context_pair_type& context_pair = contexts[tails[antecedent_index]];
	    
	    buffer.insert(buffer.end(), context_pair.first.begin(), context_pair.first.end());
	    
	    if (! context_pair.second.empty()) {
	      star_last = buffer.size();
	      if (star_first < 0)
		star_first = buffer.size();
	      
	      buffer.insert(buffer.end(), context_pair.second.begin(), context_pair.second.end());
	    }
	    
	  } else if (*citer != vocab_type::EPSILON)
	    buffer.push_back(*citer);
	}
	
	state_type state;
	if (star_first >= 0) {
	  const int prefix_size = utils::bithack::min(star_first, context_size);
	  const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
	  
	  state = std::make_pair(context_type(buffer.begin(), buffer.begin() + prefix_size),
				 context_type(buffer.end() - suffix_size, buffer.end()));
	  
	} else {
	  state = (static_cast<int>(buffer.size()) <= context_size
		   ? std::make_pair(context_type(buffer.begin(), buffer.end()),
				    context_type())
		   : std::make_pair(context_type(buffer.begin(), buffer.begin() + context_size),
				    context_type(buffer.end() - context_size, buffer.end())));
	}
	
	return state;
      }
    }

    
  private:
    context_pair_set_type contexts;
    
    node_map_type  node_map;
    state_set_type state_buf;
    
    buffer_type buffer;
    
    int  order;
  };
  
  namespace impl
  {
    template <typename Counts>
    struct expand_ngram_op
    {
      template <typename Edge, typename Weight>
      void operator()(const Edge& edge, const Weight& weight, Counts& counts) const
      {
	// no op
      }

      template <typename Edge, typename Weight, typename Iterator>
      void operator()(const Edge& edge, const Weight& weight, Counts& counts, Iterator first, Iterator last) const
      {
	counts[typename Counts::key_type(first, last)] += weight;
      }
    };
  };
  
  
  inline
  void expand_ngram(const HyperGraph& source, HyperGraph& target, const int order)
  {
    ExpandNGram __expand(order);
    
    __expand(source, target);
  }

  
};

#endif
