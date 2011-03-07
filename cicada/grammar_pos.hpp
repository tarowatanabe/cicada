// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__GRAMMAR_POS__HPP__
#define __CICADA__GRAMMAR_POS__HPP__ 1

// very siple mutable grammar class..

#include <string>
#include <vector>

#include <cicada/grammar_mutable.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/lattice.hpp>

#include <google/dense_hash_set>


namespace cicada
{
  class GrammarPOS : public GrammarPOS
  {
  public:
    typedef Lattice    lattice_type;
    
  public:
    GrammarPOS(const lattice_type& lattice)
      : GrammarMutable(1)
    {
      typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > symbol_set_type;
      
      symbol_set_type symbols;
      symbols.set_empty_key(symbol_type());
      
      feature_set_type features;
      features["pos"] = - 1.0;
      
      for (size_t first = 0; first != lattice.size(); ++ first) {
	const lattice_type::arc_set_type& arcs = lattice[first];
	
	lattice_type::arc_set_type::const_iterator aiter_end = arcs.end();
	for (lattice_type::arc_set_type::const_iterator aiter = arcs.begin(); aiter != aiter_end; ++ aiter)
	  if (aiter->label != vocab_type::EPSILON && symbols.insert(aiter->label).second) {
	    const symbol_type terminal = aiter->label.terminal();
	    const symbol_type pos = aiter->label.pos();
	    
	    if (pos == symbol_type())
	      throw std::runtime_error("no pos? " + static_cast<const std::string&>(aiter->label));
	    
	    rule_ptr_type rule(rule_type::create(rule_type(pos, rule_type::symbol_set_type(1, terminal))));
	    
	    insert(rule, rule, features, attributes);
	  }
      }
    }
  };
};

#endif
