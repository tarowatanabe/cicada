// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__REMOVE_EPSILON__HPP__
#define __CICADA__REMOVE_EPSILON__HPP__ 1

#include <cicada/vocab.hpp>
#include <cicada/remove_symbol.hpp>
  
namespace cicada
{
  namespace detail
  {
    typedef cicada::Vocab vocab_type;
    typedef vocab_type::symbol_type symbol_type;
    
    struct remove_epsilon
    {
      bool operator()(const symbol_type& x) const
      {
	return x == vocab_type::EPSILON;
      }
    };
  };

  inline
  void remove_epsilon(Lattice& lattice)
  {
    RemoveSymbol __remover;
    Lattice __lattice;
    __remover(lattice, __lattice, detail::remove_epsilon());
    lattice.swap(__lattice);
  }
  
  inline
  void remove_epsilon(const Lattice& lattice, Lattice& removed)
  {
    RemoveSymbol __remover;
    __remover(lattice, removed, detail::remove_epsilon());
  }

  inline
  void remove_epsilon(HyperGraph& graph)
  {
    RemoveSymbol __remover;
    HyperGraph removed;
    __remover(graph, removed, detail::remove_epsilon());
    graph.swap(removed);
  }
  
  inline
  void remove_epsilon(const HyperGraph& graph, HyperGraph& removed)
  {
    RemoveSymbol __remover;
    __remover(graph, removed, detail::remove_epsilon());
  }
  
};


#endif
