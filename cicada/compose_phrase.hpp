// -*- mode: c++ -*-

#ifndef __CICADA__COMPOSE_PHRASE__HPP__
#define __CICADA__COMPOSE_PHRASE__HPP__ 1

#include <vector>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/hashmurmur.hpp>

#include <google/dense_hash_map>

namespace cicada
{
  struct ComposePhrase
  {
    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;

    typedef Lattice    lattice_type;
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    ComposePhrase(const grammar_type& __grammar)
      : grammar(__grammar)
    {
      // initializer...
      
      rule_epsilon.reset(new rule_type(vocab_type::X,
				       rule_type::symbol_set_type(1, vocab_type::EPSILON),
				       rule_type::symbol_set_type(1, vocab_type::EPSILON)));
      
      rule_goal.reset(new rule_type(vocab_type::GOAL,
				    rule_type::symbol_set_type(1, vocab_type::X1),
				    rule_type::symbol_set_type(1, vocab_type::X1),
				    1));
      
      rule_x1.reset(new rule_type(vocab_type::X,
				  rule_type::symbol_set_type(1, vocab_type::X1),
				  rule_type::symbol_set_type(1, vocab_type::X1),
				  1));
      
      std::vector<symbol_type, std::allocator<symbol_type> > sequence(2);
      sequence.front() = vocab_type::X1;
      sequence.back()  = vocab_type::X2;
      
      rule_x1_x2.reset(new rule_type(vocab_type::X,
				     rule_type::symbol_set_type(sequence.begin(), sequence.end()),
				     rule_type::symbol_set_type(sequence.begin(), sequence.end()),
				     2));
    }

    void operator()(const lattice_type& lattice, hypergraph_type& graph)
    {
      
    }
    
    
    const grammar_type& grammar;

    rule_ptr_type rule_goal;
    rule_ptr_type rule_epsilon;
    rule_ptr_type rule_x1;
    rule_ptr_type rule_x1_x2;
    
  };
  
  inline
  void compose_phrase(const Grammar& grammar, const Lattice& lattice, HyperGraph& graph)
  {
    ComposePhrase __composer(grammar);
    __composer(lattice, graph);
  }

};

#endif
