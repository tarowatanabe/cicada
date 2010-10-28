// -*- mode: c++ -*-

#ifndef __CICADA__TREE_RULE__HPP__
#define __CICADA__TREE_RULE__HPP__ 1

#include <iostream>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/symbol_vector.hpp>
#include <cicada/feature_vector.hpp>

namespace cicada
{
  class TreeRule
  {
  public:
    typedef cicada::Symbol       symbol_type;
    typedef cicada::Symbol       label_type;
    typedef cicada::Vocab        vocab_type;
    
    typedef utils::simple_vector<TreeRule, std::allocator<TreeRule> > antecedent_set_type;

    typedef antecedent_set_type::const_iterator const_iterator;

  public:
    TreeRule() : label(), antecedents() {}
    explicit TreeRule(const label_type& __label) : label(__label), antecedents() {}
    explicit TreeRule(const std::string& x) : label(), antecedents() { assign(x); }
    explicit TreeRule(const char* x) : label(), antecedents() { assign(x); }
    
  public:
    void assign(const std::string& x);
    bool assign(std::string::const_iterator& iter, std::string::const_iterator end);
    
    void clear()
    {
      label = label_type();
      antecedents.clear();
    }
    
    
  public:
    const_iterator begin() const { return antecedents.begin(); }
    const_iterator end()   const { return antecedents.end(); }

  public:
    friend
    std::istream& operator>>(std::istream& is, TreeRule& x);
    friend
    std::ostream& operator<<(std::ostream& os, const TreeRule& x);
    
  public:
    label_type         label;
    antecedent_set_type antecedents;
  };
};


#endif
