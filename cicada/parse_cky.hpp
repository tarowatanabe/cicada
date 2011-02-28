// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_CKY__HPP__
#define __CICADA__PARSE_CKY__HPP__ 1

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
  // semiring and function to compute semiring from a feature vector
  // CKY + cube-pruning (Algorithm 1)
  // 

  template <typename Semiring, typename Function>
  struct ParseCKY
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Lattice    lattice_type;
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    typedef HyperGraph hypergraph_type;

    typedef Semiring semiring_type;
    typedef Semiring score_type;
    
    typedef Function function_type;
    
    typedef hypergraph_type::id_type   id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;
    
    typedef attribute_set_type::attribute_type attribute_type;

    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;
    
    ParseCKY(const symbol_type& __goal,
	     const grammar_type& __grammar,
	     const function_type& __function,
	     const int __cube_size_max,
	     const bool __yield_source=false,
	     const bool __treebank=false)
      : goal(__goal), grammar(__grammar), function(__function), cube_size_max(__cube_size_max), yield_source(__yield_source), treebank(__treebank),
	attr_span_first("span-first"),
	attr_span_last("span-last")
    {
      goal_rule = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, goal.non_terminal())));
      
      node_map.set_empty_key(symbol_level_type());
    }
    
    struct ActiveItem
    {
      ActiveItem(const hypergraph_type::id_type& __node,
		 const hypergraph_type::edge_type::node_set_type __tails,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node),
	  tails(__tails),
	  features(__features),
	  attributes(__attributes) {}
      ActiveItem(const hypergraph_type::id_type& __node,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node),
	  tails(),
	  features(__features),
	  attributes(__attributes) {}
      ActiveItem(const hypergraph_type::id_type& __node)
	: node(__node),
	  tails(),
	  features(),
	  attributes() {}
      
      hypergraph_type::id_type                  node;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type                          features;
      attribute_set_type                        attributes; 
      score_type                                score;
    };
    
    typedef ActiveItem active_type;
    typedef utils::chunk_vector<active_type, 4096 / sizeof(active_type), std::allocator<active_type> > active_set_type;
    
    // we need candidate structure, a pointer to active-item, actual rule and 
  };
  
};

#endif
