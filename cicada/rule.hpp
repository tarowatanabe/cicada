// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__RULE__HPP__
#define __CICADA__RULE__HPP__ 1

#include <iostream>
#include <vector>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/symbol_vector.hpp>

#include <utils/bithack.hpp>
#include <utils/piece.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  class Rule
  {
  public:
    typedef cicada::Symbol       symbol_type;
    typedef cicada::Symbol       word_type;
    typedef cicada::Vocab        vocab_type;
    typedef cicada::SymbolVector symbol_set_type;

    typedef boost::shared_ptr<Rule> rule_ptr_type;
    
  public:
    Rule() : lhs(), rhs() {}
    Rule(const utils::piece& x) { assign(x); }
    Rule(const symbol_type& x_lhs, const symbol_set_type& x_rhs) : lhs(x_lhs), rhs(x_rhs) {}
    template <typename Iterator>
    Rule(const symbol_type& __lhs, Iterator first, Iterator last) : lhs(__lhs), rhs(first, last) {}

  public:
    static rule_ptr_type create(const Rule& x);
    
  public:
    void assign(const utils::piece& x);
    bool assign(std::string::const_iterator& iter, std::string::const_iterator end);
    
    void clear()
    {
      lhs = symbol_type();
      rhs.clear();
    }
    
    friend
    std::ostream& operator<<(std::ostream& os, const Rule& x);
    friend
    std::istream& operator>>(std::istream& is, Rule& x);
    
  public:
    symbol_type      lhs;
    symbol_set_type  rhs;
  };
  
  inline
  size_t hash_value(Rule const& x)
  {
    return utils::hashmurmur<size_t>()(x.rhs.begin(), x.rhs.end(), x.lhs.id());
  }

  inline
  void sort(Rule& x, Rule& y)
  {
    // sort, so that x's non-terminals are indexed in sorted order...
    typedef std::vector<int, std::allocator<int> > index_type;
    
    Rule::symbol_set_type x_new(x.rhs);
    Rule::symbol_set_type y_new(y.rhs);
    index_type index(utils::bithack::max(x.rhs.size(), y.rhs.size()) + 1);
    
    int pos = 1;
    Rule::symbol_set_type::iterator xiter_end = x_new.end();
    for (Rule::symbol_set_type::iterator xiter = x_new.begin(); xiter != xiter_end; ++ xiter) 
      if (xiter->is_non_terminal()) {
	const int non_terminal_pos = xiter->non_terminal_index();
	
	index[utils::bithack::branch(non_terminal_pos == 0, pos, non_terminal_pos)] = pos;
	
	*xiter = xiter->non_terminal(pos);
	++ pos;
      }
    
    pos = 1;
    Rule::symbol_set_type::iterator yiter_end = y_new.end();
    for (Rule::symbol_set_type::iterator yiter = y_new.begin(); yiter != yiter_end; ++ yiter)
      if (yiter->is_non_terminal()) {
	*yiter = yiter->non_terminal(index[pos]);
	++ pos;
      }
    
    x.rhs.swap(x_new);
    y.rhs.swap(y_new);
  }


  inline
  bool operator==(const Rule& x, const Rule& y)
  {
    return x.lhs == y.lhs && x.rhs == y.rhs;
  }

  inline
  bool operator!=(const Rule& x, const Rule& y)
  {
    return x.lhs != y.lhs || x.rhs != y.rhs;
  }

  inline
  bool operator<(const Rule& x, const Rule& y)
  {
    return x.lhs < y.lhs || (!(y.lhs < x.lhs) && x.rhs < y.rhs);
  }
  
  inline
  bool operator>(const Rule& x, const Rule& y)
  {
    return y < x;
  }

  inline
  bool operator<=(const Rule& x, const Rule& y)
  {
    return ! (y < x);
  }
  
  inline
  bool operator>=(const Rule& x, const Rule& y)
  {
    return ! (x < y);
  }
  
};

#endif
