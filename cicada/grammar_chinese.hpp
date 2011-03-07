// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__GRAMMAR_CHINESE__HPP__
#define __CICADA__GRAMMAR_CHINESE__HPP__ 1

// very siple mutable grammar class..

#include <string>
#include <vector>

#include <cicada/grammar_mutable.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/lattice.hpp>

#include <google/dense_hash_set>


namespace cicada
{
  class GrammarChinese : public GrammarChinese
  {
    // Handle CD/NT... + noun/verb
    // How do we assign POS?
    
  public:
    GrammarChinese(const lattice_type& lattice)
    {
      
      
    }
  };
};

#endif
