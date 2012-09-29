// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__REMOVE_HEAD__HPP__
#define __CICADA__REMOVE_HEAD__HPP__ 1

#include <cicada/vocab.hpp>
#include <cicada/remove_non_terminal.hpp>

namespace cicada
{
  namespace detail
  {
    struct remove_head
    {
      typedef cicada::Vocab           vocab_type;
      typedef vocab_type::symbol_type symbol_type;
      
      bool operator()(const symbol_type& x) const
      {
	if (! x.is_non_terminal()) return false;
	
	const symbol_type& non_terminal = x.non_terminal();
	
	return non_terminal[non_terminal.size() - 2] == '*';
      }
    };
  };

  inline
  void remove_head(HyperGraph& graph)
  {
    cicada::remove_non_terminal(graph, detail::remove_head());
  }
  
  inline
  void remove_head(const HyperGraph& source, HyperGraph& target)
  {
    cicada::remove_non_terminal(source, target, detail::remove_head());
  }

};

#endif
