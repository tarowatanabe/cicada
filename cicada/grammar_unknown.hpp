// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__GRAMMAR_UNKNOWN__HPP__
#define __CICADA__GRAMMAR_UNKNOWN__HPP__ 1

#include <string>
#include <vector>

#include <cicada/grammar.hpp>
#include <cicada/grammar_mutable.hpp>
#include <cicada/signature.hpp>

#include <utils/compact_map.hpp>
#include <utils/compact_set.hpp>
#include <utils/trie_compact.hpp>
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

  public:
    typedef Grammar::transducer_ptr_type transducer_ptr_type;

  private:

    struct unassigned_unigram
    {
      uchar_type operator()() const { return 0; }
    };


    typedef utils::trie_compact<symbol_type, double,
				utils::unassigned<symbol_type>, 
				boost::hash<symbol_type>, std::equal_to<symbol_type>,
				std::allocator<std::pair<const symbol_type, double> > > backoff_set_type;
    
    typedef utils::compact_map<uchar_type, double,
			       unassigned_unigram, unassigned_unigram,
			       utils::hashmurmur<size_t>, std::equal_to<uchar_type>,
			       std::allocator<std::pair<const uchar_type, double> > > unigram_set_type;
    
    typedef utils::trie_compact<symbol_type, unigram_set_type,
				utils::unassigned<symbol_type>, 
				boost::hash<symbol_type>, std::equal_to<symbol_type>,
				std::allocator<std::pair<const symbol_type, unigram_set_type> > > ngram_set_type;
    
  public:
    GrammarUnknown(const std::string& __signature,
		   const std::string& __parameter)
      : base_type(1),
	signature(&signature_type::create(__signature)),
	__grammar_oov(new GrammarMutable(1)),
	backoff(),
	ngram(),
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
	__grammar_oov(new GrammarMutable(1)),
	backoff(),
	ngram(),
	unigram(),
	logprob_unk(0),
	feature_character("rule-character")
    {
      base_type::read(__parameter);
      read_character(__character);
    }
    
    transducer_ptr_type clone() const
    {
      std::auto_ptr<GrammarUnknown> __tmp(new GrammarUnknown(*this));
      
      __tmp->signature = &signature_type::create(signature->algorithm());
      if (__grammar_oov)
	__tmp->__grammar_oov = __grammar_oov->clone();
      
      return transducer_ptr_type(__tmp.release());
    }

    transducer_ptr_type grammar_oov() const
    {
      return __grammar_oov;
    }
    
    void assign(const hypergraph_type& graph)
    {
      GrammarMutable* oov = dynamic_cast<GrammarMutable*>(&(*__grammar_oov));
      
      if (! oov)
	throw std::runtime_error("no oov grammar?");

      if (oov->size() > 1024 * 4)
	oov->clear();

      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) 
	if (eiter->rule) {
	  const rule_type& rule = *(eiter->rule);
	  
	  rule_type::symbol_set_type::const_iterator siter_end = rule.rhs.end();
	  for (rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter) 
	    if (*siter != vocab_type::EPSILON && siter->is_terminal())
	      if (base_type::next(base_type::root(), *siter) == base_type::root())
		insert(*siter);
	}
    }
    
    void assign(const lattice_type& lattice)
    {
      GrammarMutable* oov = dynamic_cast<GrammarMutable*>(&(*__grammar_oov));
      
      if (! oov)
	throw std::runtime_error("no oov grammar?");
      
      if (oov->size() > 1024 * 4)
	oov->clear();

      for (size_t first = 0; first != lattice.size(); ++ first) {
	const lattice_type::arc_set_type& arcs = lattice[first];
	
	lattice_type::arc_set_type::const_iterator aiter_end = arcs.end();
	for (lattice_type::arc_set_type::const_iterator aiter = arcs.begin(); aiter != aiter_end; ++ aiter)
	  if (aiter->label != vocab_type::EPSILON)
	    if (base_type::next(base_type::root(), aiter->label) == base_type::root())
	      insert(aiter->label);
      }
    }
    
  private:
    void insert(const symbol_type& word);
    void read_character(const std::string& file);
    
  private:
    const signature_type* signature;
    
    transducer_ptr_type __grammar_oov;

    backoff_set_type backoff;
    ngram_set_type   ngram;
    unigram_set_type unigram;
    
    double logprob_unk;
    
    feature_type feature_character;
  };
};

#endif
