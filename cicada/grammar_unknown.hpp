// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__GRAMMAR_UNKNOWN__HPP__
#define __CICADA__GRAMMAR_UNKNOWN__HPP__ 1

#include <string>
#include <vector>

#include <cicada/grammar_mutable.hpp>
#include <cicada/signature.hpp>

#include <google/dense_hash_map>

#include <utils/compact_trie_dense.hpp>
#include <utils/hashmurmur.hpp>

namespace cicada
{
  class GrammarUnknown : public GrammarMutable
  {
  private:
    typedef GrammarMutable base_type;
    typedef Signature signature_type;
    
    typedef int32_t uchar_type;

    typedef feature_set_type::feature_type feature_type;

  private:
    typedef utils::compact_trie_dense<symbol_type, double, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				      std::allocator<std::pair<const symbol_type, double> > > backoff_set_type;
    class unigram_set_type : public google::dense_hash_map<uchar_type, double, utils::hashmurmur<size_t>, std::equal_to<uchar_type> >
    {
    public:
      typedef google::dense_hash_map<uchar_type, double, utils::hashmurmur<size_t>, std::equal_to<uchar_type> > base_type;
      
      unigram_set_type() : base_type() { base_type::set_empty_key(0); }
    };
    
    typedef utils::compact_trie_dense<symbol_type, unigram_set_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				      std::allocator<std::pair<const symbol_type, unigram_set_type> > > ngram_set_type;
    
    
  public:
    GrammarUnknown(const std::string& __signature,
		   const std::string& __parameter)
      : base_type(1),
	signature(&signature_type::create(__signature)),
	backoff(symbol_type()),
	ngram(symbol_type()),
	unigram(),
	logprob_unk(0),
	feature_character()
    {
      base_type::read(__parameter);
    }
    
    GrammarUnknown(const std::string& __signature,
		   const std::string& __parameter,
		   const std::string& __character)
      : base_type(1),
	signature(&signature_type::create(__signature)),
	backoff(symbol_type()),
	ngram(symbol_type()),
	unigram(),
	logprob_unk(0),
	feature_character("character-penalty")
    {
      base_type::read(__parameter);
      read_character(__character);
    }
    
    transducer_ptr_type clone() const
    {
      std::auto_ptr<GrammarUnknown> __tmp(new GrammarUnknown(*this));
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
    void insert(const symbol_type& word);
    void read_character(const std::string& file);
    
  private:
    const signature_type* signature;
    
    backoff_set_type backoff;
    ngram_set_type   ngram;
    unigram_set_type unigram;
    
    double logprob_unk;
    
    feature_type feature_character;
  };
};

#endif
