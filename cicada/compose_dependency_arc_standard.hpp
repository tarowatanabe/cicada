// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__COMPOSE_DEPENDENCY_ARC_STANDARD__HPP__
#define __CICADA__COMPOSE_DEPENDENCY_ARC_STANDARD__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/bit_vector.hpp>

#include <utils/sgi_hash_set.hpp>
#include <utils/sgi_hash_map.hpp>

namespace cicada
{
  // arc standard...!
  struct ComposeDependencyArcStandard
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

    ComposeDependencyArcStandard(const grammar_tyep& __grammar,
				 const bool __pos_mode=false)
      : grammar(__grammar), pos_mode(__pos_mode),
	attr_dependency_head("dependency-head"),
	attr_dependency_dependent("dependency-dependent") {}
    
    void operator()(const lattice_type& lattice,
		    hypegraph_type& graph)
    {
      
      
    }
    
  private:
    const gramamr_type& grammar;
    const bool pos_mode;

    const attribute_type attr_dependency_head;
    const attribute_type attr_dependency_dependent;
    
    
  };
};

#endif
