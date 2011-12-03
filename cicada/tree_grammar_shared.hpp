// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TREE_GRAMMAR_SHARED__HPP__
#define __CICADA__TREE_GRAMMAR_SHARED__HPP__ 1

// shared storage grammar...

#include <string>

#include <cicada/tree_transducer.hpp>
#include <cicada/tree_grammar_mutable.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  
  class TreeGrammarShared : public TreeTransducer
  {
  private:
    typedef TreeGrammarMutable impl_type;
    
  public:
    TreeGrammarShared(const std::string& parameter) : pimpl(new impl_type(parameter)) {}
    ~TreeGrammarShared() {}
    
    TreeGrammarShared(const TreeGrammarShared& x) : pimpl(x.pimpl) {}
    TreeGrammarShared& operator=(const TreeGrammarShared& x)
    {
      pimpl = x.pimpl;
      return *this;
    }

  private:
    TreeGrammarShared() {}
    
  public:
    // virtual members
    transducer_ptr_type clone() const { return transducer_ptr_type(new TreeGrammarShared(*this)); }
    
    edge_type edge(const symbol_type& symbol) const { return pimpl->edge(symbol); }
    edge_type edge(const symbol_set_type& symbols) const { return pimpl->edge(symbols); }
    edge_type edge(const symbol_type* first, const symbol_type* last) const { return pimpl->edge(first, last); }
    
    id_type root() const { return pimpl->root(); }
    id_type next(const id_type& node, const edge_type& edge) const { return pimpl->next(node, edge); }
    id_type next(const id_type& node, const symbol_type& symbol) const { return pimpl->next(node, symbol); }
    bool has_next(const id_type& node) const { return pimpl->has_next(node); }
    const rule_pair_set_type& rules(const id_type& node) const { return pimpl->rules(node); }

  private:
    boost::shared_ptr<impl_type> pimpl;
  };
  
};

#endif
