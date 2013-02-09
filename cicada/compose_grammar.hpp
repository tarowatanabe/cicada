// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_GRAMMAR__HPP__
#define __CICADA__COMPOSE_GRAMMAR__HPP__ 1

#include <vector>
#include <algorithm>
#include <set>
#include <sstream>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/unordered_map.hpp>

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
    
    struct rule_hash_type
    {
      size_t operator()(const rule_type* x) const
      {
	return (x ? hash_value(*x) : size_t(0));
      }
    };
    
    struct rule_equal_type
    {
      bool operator()(const rule_type* x, const rule_type* y) const
      {
	return x == y ||(x && y && *x == *y);
      }
    };

    typedef utils::unordered_map<const rule_type*, std::string, rule_hash_type, rule_equal_type,
				 std::allocator<std::pair<const rule_type*, std::string> > >::type frontier_set_type;
    
    ComposeGrammar(const grammar_type& __grammar, const bool __yield_source=false, const bool __frontier=false)
      : grammar(__grammar),
	yield_source(__yield_source),
	frontier(__frontier),
	attr_frontier_source(__frontier ? "frontier-source" : ""),
        attr_frontier_target(__frontier ? "frontier-target" : "")
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

      frontiers_source.clear();
      frontiers_target.clear();
      
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

	    if (frontier) {
	      const rule_type* rule_source = riter->source.get();
	      const rule_type* rule_target = riter->target.get();
	      
	      if (rule_source) {
		frontier_set_type::iterator siter = frontiers_source.find(rule_source);
		if (siter == frontiers_source.end()) {
		  std::ostringstream os;
		  os << rule_source->rhs;
		  
		  siter = frontiers_source.insert(std::make_pair(rule_source, os.str())).first;
		}
		
		edge_new.attributes[attr_frontier_source] = siter->second;
	      }
	      
	      if (rule_target) {
		frontier_set_type::iterator titer = frontiers_target.find(rule_target);
		if (titer == frontiers_target.end()) {
		  std::ostringstream os;
		  os << rule_target->rhs;
		  
		  titer = frontiers_target.insert(std::make_pair(rule_target, os.str())).first;
		}
		
		edge_new.attributes[attr_frontier_target] = titer->second;
	      }
	    }
	    
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
    const bool frontier;
    
    const attribute_type attr_frontier_source;
    const attribute_type attr_frontier_target;

    frontier_set_type frontiers_source;
    frontier_set_type frontiers_target;
  };

  inline
  void compose_grammar(const Grammar& grammar, const HyperGraph& source, HyperGraph& target, const bool yield_source=false, const bool frontier=false)
  {
    ComposeGrammar(grammar, yield_source, frontier)(source, target);
  }

};

#endif
