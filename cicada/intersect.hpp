// -*- mode: c++ -*-

#ifndef __CICADA__INTERSECT__HPP__
#define __CICADA__INTERSECT__HPP__ 1

#include <stdexcept>

#include <cicada/rule.hpp>
#include <cicada/compose_cky.hpp>
#include <cicada/grammar.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar_mutable.hpp>
#include <cicada/vocab.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  inline
  void intersect(const HyperGraph& source, const Lattice& lattice, HyperGraph& target)
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Vocab      vocab_type;
    typedef Symbol     symbol_type;
    typedef HyperGraph hypergraph_type;
    typedef Lattice    lattice_type;
    
    typedef Grammar        grammar_type;
    typedef GrammarMutable grammar_mutable_type;

    typedef hypergraph_type::node_type     node_type;
    typedef hypergraph_type::edge_type     edge_type;
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef grammar_type::rule_pair_type rule_pair_type;
    
    typedef rule_type::symbol_set_type  symbol_set_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > yield_type;

    target.clear();
    
    if (source.goal == hypergraph_type::invalid || lattice.empty()) return;
    
    non_terminal_set_type non_terminals(source.nodes.size());
    for (size_type id = 0; id < source.nodes.size(); ++ id)
      non_terminals[id] = std::string("[NODE_") + boost::lexical_cast<std::string>(id) + ']';
    
    
    boost::shared_ptr<grammar_mutable_type> g(new grammar_mutable_type());
    grammar_type grammar;
    grammar.push_back(g);

    yield_type yield;
    
    hypergraph_type::edge_set_type::const_iterator eiter_end = source.edges.end();
    for (hypergraph_type::edge_set_type::const_iterator eiter = source.edges.begin(); eiter != eiter_end; ++ eiter) {
      const edge_type& edge = *eiter;
      
      yield.clear();
      yield.insert(yield.end(), edge.rule->rhs.begin(), edge.rule->rhs.end());
      
      
      // sort...
      if (! edge.tails.empty()) {
	int pos = 0;
	yield_type::iterator siter_end = yield.end();
	for (yield_type::iterator siter = yield.begin(); siter != siter_end; ++ siter)
	  if (siter->is_non_terminal()) {
	    int non_terminal_index = siter->non_terminal_index() - 1;
	    if (non_terminal_index < 0)
	      non_terminal_index = pos;
	    
	    *siter = non_terminals[edge.tails[non_terminal_index]].non_terminal(pos);
	    ++ pos;
	  }
      }
      
      rule_ptr_type rule(new rule_type(non_terminals[edge.head], symbol_set_type(yield.begin(), yield.end())));
      
      g->insert(rule_pair_type(rule, rule, edge.features));
    }
    
    const symbol_type goal = non_terminals[source.goal];
    
    ComposeCKY composer(goal, grammar);
    
    composer(lattice, target);
  }
};

#endif
