// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__REMOVE_BOS_EOS__HPP__
#define __CICADA__REMOVE_BOS_EOS__HPP__ 1

#include <cicada/vocab.hpp>
#include <cicada/remove_symbol.hpp>

namespace cicada
{

  namespace detail
  {
    typedef cicada::Vocab vocab_type;
    typedef vocab_type::symbol_type symbol_type;
    
    struct remove_bos_eos
    {
      bool operator()(const symbol_type& x) const
      {
	return x == vocab_type::BOS || x == vocab_type::EOS || x == vocab_type::EPSILON;
      }
    };
  };
  
  inline
  void remove_bos_eos(Lattice& lattice)
  {
    RemoveSymbol __remover;
    Lattice __lattice;
    __remover(lattice, __lattice, detail::remove_bos_eos());
    lattice.swap(__lattice);
  }
  
  inline
  void remove_bos_eos(const Lattice& lattice, Lattice& removed)
  {
    RemoveSymbol __remover;
    __remover(lattice, removed, detail::remove_bos_eos());
  }
  
  inline
  void remove_bos_eos(HyperGraph& graph)
  {
    RemoveSymbol __remover;
    HyperGraph removed;
    __remover(graph, removed, detail::remove_bos_eos());
    graph.swap(removed);
  }
  
  inline
  void remove_bos_eos(const HyperGraph& graph, HyperGraph& removed)
  {
    RemoveSymbol __remover;
    __remover(graph, removed, detail::remove_bos_eos());
  }
  
};


#endif
