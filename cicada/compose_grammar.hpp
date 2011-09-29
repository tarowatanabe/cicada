// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_GRAMMAR__HPP__
#define __CICADA__COMPOSE_GRAMMAR__HPP__ 1

#include <vector>
#include <algorithm>
#include <set>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/sgi_hash_map.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

namespace cicada
{
  //
  // compose hypergraph with grammar, yielding another hypergraph!
  //
  struct ComposeGrammar
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;

    typedef Lattice    lattice_type;
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    typedef HyperGraph hypergraph_type;
  
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef attribute_set_type::attribute_type attribute_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;
    
    
    ComposeGrammar(const grammar_type& __grammar, const bool __yield_source=false)
      : grammar(__grammar), yield_source(__yield_source),
    { }

    struct filter_edge
    {
      typedef std::vector<bool, std::allocator<bool> > removed_type;
      
      filter_edge(const removed_type& __removed) : removed(__removed) {}
  
      bool operator()(const hypergraph_type::edge_type& edge) const
      {
	return removed[edge.id];
      }
      
      const removed_type& removed;
    };
    
    void operator()(const hypergraph_type& source, hypergraph_type& graph)
    {
      graph = source;
      if (! graph.is_valid()) return;

      filter_edge::removed_type removed(graph.edges.size(), false);
      bool found = false;
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = *eiter;
	
	for (size_type table = 0; table != grammar.size(); ++ table) {
	  const transducer_type& transducer = grammar[table];

	  transducer_type::id_type node = transducer.root();
	  
	  {
	    rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	    for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) {
	      node = transducer.next(node, *riter);
	      if (node == transducer.root()) break;
	    }
	  }
	  
	  if (node == transducer.root()) continue;
	  
	  const transducer_type::rule_pair_set_type& rules = transducer.rules(node);
	  
	  if (rules.empty()) continue;
	  
	  transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	  for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	    const rule_ptr_type rule = (yield_source ? riter->source : riter->target);
	    const symbol_type& lhs = rule->lhs;
	    
	    if (lhs != edge.rule->lhs) continue;
	    
	    // if we found mathcing, simply replace with the new one!
	    removed[edge.id] = true;
	    found = true;
	    
	    // TODO: we need to sort source-side, but we will not check here...
	    hypergraph_type::edge_type& edge_new = graph.add_edge(edge.tails.begin(), edge.tails.end());
	    edge_new.rule = rule;
	    edge_new.features = riter->features;
	    edge_new.attributes = riter->attributes;
	    
	    graph.connect_edge(edge_new.id, edge.head);
	  }
	}
      }
      
      removed.resize(graph.edges.size());
      
      if (found) {
	hypergraph_type sorted;
	topologically_sort(graph, sorted, filter_edge(removed), true);
	graph.swap(sorted);
      }
    }
    
  private:
    const grammar_type& grammar;
    const bool yield_source;
  };

  inline
  void compose_grammar(const Grammar& grammar, const HyperGraph& source, HyperGraph& target, const bool yield_source=false)
  {
    ComposeGrammar(grammar, yield_source)(source, target);
  }

};

#endif
