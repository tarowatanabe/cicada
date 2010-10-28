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

    typedef antecedent_set_type::size_type       size_type;
    typedef antecedent_set_type::difference_type difference_type;
    
    typedef antecedent_set_type::const_iterator  const_iterator;
    typedef antecedent_set_type::iterator        iterator;
    typedef antecedent_set_type::const_reference const_reference;
    typedef antecedent_set_type::reference       reference;
    
  public:
    TreeRule() : label(), antecedents() {}
    TreeRule(const label_type& __label) : label(__label), antecedents() {}
    template <typename Iterator>
    TreeRule(const label_type& __label, Iterator first, Iterator last) : label(__label), antecedents(first, last) {}
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
    iterator begin() { return antecedents.begin(); }
    
    const_iterator end()   const { return antecedents.end(); }
    iterator end() { return antecedents.end(); }
    
    const_reference operator[](size_type pos) const { return antecedents[pos]; }
    reference operator[](size_type pos) { return antecedents[pos]; }
    
    const_reference front() const { return antecedents.front(); }
    reference front() { return antecedents.front(); }
    
    const_reference back() const { return antecedents.back(); }
    reference back() { return antecedents.back(); }
    

    size_type size() const { return antecedents.size(); }
    bool empty() const { return antecedents.empty(); }

    void swap(TreeRule& x)
    {
      label.swap(x.label);
      antecedents.swap(x.antecedents);
    }

  public:
    friend
    std::istream& operator>>(std::istream& is, TreeRule& x);
    friend
    std::ostream& operator<<(std::ostream& os, const TreeRule& x);
    
  public:
    label_type          label;
    antecedent_set_type antecedents;
  };
  
  inline
  size_t hash_value(TreeRule const& x)
  {
    size_t seed = utils::hashmurmur<size_t>()(x.label.id());
    for (TreeRule::const_iterator aiter = x.begin(); aiter != x.end(); ++ aiter)
      seed = utils::hashmurmur<size_t>()(hash_value(*aiter), seed);
    return seed;
  }

  inline
  bool operator==(const TreeRule& x, const TreeRule& y)
  {
    return x.label == y.label && x.antecedents == y.antecedents;
  }

  inline
  bool operator!=(const TreeRule& x, const TreeRule& y)
  {
    return x.label != y.label || x.antecedents != y.antecedents;
  }

  inline
  bool operator<(const TreeRule& x, const TreeRule& y)
  {
    return x.label < y.label || (!(y.label < x.label) && x.antecedents < y.antecedents);
  }
  
  inline
  bool operator>(const TreeRule& x, const TreeRule& y)
  {
    return y < x;
  }
  
  inline
  bool operator<=(const TreeRule& x, const TreeRule& y)
  {
    return ! (y < x);
  }
  
  inline
  bool operator>=(const TreeRule& x, const TreeRule& y)
  {
    return ! (x < y);
  }
  
};

namespace std
{
  inline
  void swap(cicada::TreeRule& x, cicada::TreeRule& y)
  {
    x.swap(y);
  }
};


#endif
