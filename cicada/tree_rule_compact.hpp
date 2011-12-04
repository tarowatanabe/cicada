// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TREE_RULE_COMPACT__HPP__
#define __CICADA__TREE_RULE_COMPACT__HPP__ 1

// a compact representation of tree-rule

#include <cicada/tree_rule.hpp>

#include <utils/simple_vector.hpp>
#include <utils/hashmurmur.hpp>

namespace cicada
{
  class TreeRuleCompact
  {
  public:
    typedef cicada::Symbol       symbol_type;
    typedef cicada::Symbol       label_type;

    typedef cicada::TreeRule tree_rule_type;
    typedef cicada::TreeRule rule_type;
    
    typedef tree_rule_type::rule_ptr_type tree_rule_ptr_type;
    typedef tree_rule_type::rule_ptr_type rule_ptr_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

  private:
    typedef uint8_t byte_type;
    typedef utils::simple_vector<byte_type, std::allocator<byte_type> > impl_type;
    
  public:
    TreeRuleCompact() {}
    TreeRuleCompact(const TreeRule& x) { assign(x); }
    TreeRuleCompact(const TreeRuleCompact& x) : impl(x.impl) {}
    
    TreeRuleCompact& operator=(const TreeRule& x)
    {
      assign(x);
      return *this;
    }
    
    TreeRuleCompact& operator=(const TreeRuleCompact& x)
    {
      impl = x.impl;
      return *this;
    }
    
    void assign(const TreeRuleCompact& x)
    {
      impl = x.impl;
    }
    
    // actual assignment...
    void assign(const TreeRule& x) { encode(x); }
    
    // encode/decode
    void encode(const TreeRule& x);
    TreeRule decode() const;
    
    size_type size_compressed() const  { return impl.size(); }
    
    void clear() { impl.clear(); }

    void swap(TreeRuleCompact& x)
    {
      impl.swap(x.impl);
    }
    
  public:
    friend size_t hash_value(TreeRuleCompact const& x) { return utils::hashmurmur<size_t>()(x.impl.begin(), x.impl.end(), 0); }
    friend bool operator==(const TreeRuleCompact& x, const TreeRuleCompact& y) { return x.impl == y.impl; }
    friend bool operator!=(const TreeRuleCompact& x, const TreeRuleCompact& y) { return x.impl != y.impl; }
    friend bool operator<(const TreeRuleCompact& x, const TreeRuleCompact& y) { return x.impl < y.impl; }
    friend bool operator>(const TreeRuleCompact& x, const TreeRuleCompact& y) { return x.impl > y.impl; }
    friend bool operator<=(const TreeRuleCompact& x, const TreeRuleCompact& y) { return x.impl <= y.impl; }
    friend bool operator>=(const TreeRuleCompact& x, const TreeRuleCompact& y) { return x.impl >= y.impl; }

  private:
    impl_type impl;
  };
};

namespace std
{
  inline
  void swap(cicada::TreeRuleCompact& x, cicada::TreeRuleCompact& y)
  {
    x.swap(y);
  }
};

#endif
