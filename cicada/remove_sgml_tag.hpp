// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__REMOVE_SGML_TAG__HPP__
#define __CICADA__REMOVE_SGML_TAG__HPP__ 1

#include <cicada/vocab.hpp>
#include <cicada/remove_symbol.hpp>

namespace cicada
{

  namespace detail
  {
    typedef cicada::Vocab vocab_type;
    typedef vocab_type::symbol_type symbol_type;
    
    struct remove_sgml_tag
    {
      bool operator()(const symbol_type& x) const
      {
	return (x != vocab_type::BOS && x != vocab_type::EOS && x.is_sgml_tag()) || x == vocab_type::EPSILON;
      }
    };

    struct remove_sgml_tag_all
    {
      bool operator()(const symbol_type& x) const
      {
	return x.is_sgml_tag() || x == vocab_type::EPSILON;
      }
    };
    
  };
  
  inline
  void remove_sgml_tag(Lattice& lattice, const bool remove_all=false)
  {
    RemoveSymbol __remover;
    Lattice __lattice;
    
    if (remove_all)
      __remover(lattice, __lattice, detail::remove_sgml_tag_all());
    else
      __remover(lattice, __lattice, detail::remove_sgml_tag());
    
    lattice.swap(__lattice);
  }
  
  inline
  void remove_sgml_tag(const Lattice& lattice, Lattice& removed, const bool remove_all=false)
  {
    RemoveSymbol __remover;
    
    if (remove_all)
      __remover(lattice, removed, detail::remove_sgml_tag_all());
    else
      __remover(lattice, removed, detail::remove_sgml_tag());
  }
  
  inline
  void remove_sgml_tag(HyperGraph& graph, const bool remove_all=false)
  {
    RemoveSymbol __remover;
    HyperGraph removed;
    
    if (remove_all)
      __remover(graph, removed, detail::remove_sgml_tag_all());
    else
      __remover(graph, removed, detail::remove_sgml_tag());
    
    graph.swap(removed);
  }
  
  inline
  void remove_sgml_tag(const HyperGraph& graph, HyperGraph& removed, const bool remove_all=false)
  {
    RemoveSymbol __remover;
    
    if (remove_all)
      __remover(graph, removed, detail::remove_sgml_tag_all());
    else
      __remover(graph, removed, detail::remove_sgml_tag());
  }
  
};


#endif
