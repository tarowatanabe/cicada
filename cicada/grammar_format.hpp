// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__GRAMMAR_FORMAT__HPP__
#define __CICADA__GRAMMAR_FORMAT__HPP__ 1

#include <string>
#include <vector>

#include <cicada/grammar_mutable.hpp>
#include <cicada/format.hpp>

#include <google/dense_hash_map>

#include <utils/compact_trie_dense.hpp>
#include <utils/hashmurmur.hpp>

namespace cicada
{
  class GrammarFormat : public GrammarMutable
  {
  private:
    typedef GrammarMutable base_type;
    typedef Signature signature_type;
    
    typedef int32_t uchar_type;
    
    typedef feature_set_type::feature_type feature_type;
    
  private:
    
  public:
    GrammarFormat(const std::string& __parameter)
      : base_type(),
	feature("format")
    {
      base_type::read(__parameter);
    }
    
    transducer_ptr_type clone() const
    {
      std::auto_ptr<GrammarFormat> __tmp(new GrammarFormat(*this));
      __tmp->signature = &signature_type::create(signature->algorithm());
      
      return transducer_ptr_type(__tmp.release());
    }
    
    void assign(const hypergraph_type& graph)
    {
      // we will collect unigrams first...
    }
    
    void assign(const lattice_type& lattice)
    {
      // we will collect unigrams first...
      
    }
    
    
  private:
    // we will keep actual rules in base_type + queried types
    
    
    feature_type feature;
  };
};

#endif
