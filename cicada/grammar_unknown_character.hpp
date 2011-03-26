// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__GRAMMAR_UNKNOWN_CHARACTER__HPP__
#define __CICADA__GRAMMAR_UNKNOWN_CHARACTER__HPP__ 1

#include <string>
#include <vector>

#include <cicada/grammar_mutable.hpp>
#include <cicada/signature.hpp>

#include <google/dense_hash_map>

namespace cicada
{
  class GrammarUnknownCharacter : public GrammarMutable
  {
  private:
    typedef GrammarMutable base_type;
    typedef Signature signature_type;
    
  public:
    GrammarUnknownCharacter(const std::string& __signature,
			    const std::string& __parameter,
			    const std::string& __character)
      : base_type(1),
	signature(&signature_type::create(__signature))
    {
      base_type::read(__parameter);
      read_character(__character);
    }
    
    transducer_ptr_type clone() const
    {
      std::auto_ptr<GrammarUnknownCharacter> __tmp(new GrammarUnknownCharacter(*this));
      __tmp->signature = &signature_type::create(signature->algorithm());
      
      return transducer_ptr_type(__tmp.release());
    }
    
    void assign(const hypergraph_type& graph)
    {
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) 
	if (eiter->rule) {
	  const rule_type& rule = *(eiter->rule);
	  
	  rule_type::symbol_set_type::const_iterator siter_end = rule.rhs.end();
	  for (rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter) 
	    if (*siter != vocab_type::EPSILON && siter->is_terminal())
	      insert(*siter);
	}
    }
    
    void assign(const lattice_type& lattice)
    {
      for (size_t first = 0; first != lattice.size(); ++ first) {
	const lattice_type::arc_set_type& arcs = lattice[first];
	
	lattice_type::arc_set_type::const_iterator aiter_end = arcs.end();
	for (lattice_type::arc_set_type::const_iterator aiter = arcs.begin(); aiter != aiter_end; ++ aiter)
	  if (aiter->label != vocab_type::EPSILON)
	    insert(aiter->label);
      }
    }
    
  private:
    void read_character(const std::string& path);
    void insert(const symbol_type& word);
    
  private:
    const signature_type* signature;
  };
};

#endif
