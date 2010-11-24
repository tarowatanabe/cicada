// -*- mode: c++ -*-

#ifndef __CICADA__TREE_RULE__HPP__
#define __CICADA__TREE_RULE__HPP__ 1

#include <iostream>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/symbol_vector.hpp>
#include <cicada/feature_vector.hpp>

#include <utils/bithack.hpp>

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

    // compute max-depth
    size_type depth() const
    {
      return __depth(*this, 0);
    }

    // convert into hyperpth
    template <typename Path>
    void hyperpath(Path& path) const
    {
      const size_type max_depth = depth();
      
      path.clear();
      path.resize(max_depth + 1);
      
      __hyperpath(*this, path, max_depth, 0);
      
      path.front().push_back(vocab_type::NONE);
    }
    
    // collect frontiers...
    template <typename Iterator>
    void frontier(Iterator iter) const
    {
      if (antecedents.empty()) {
	*iter = label;
	++ iter;
      }
      
      const_iterator aiter_end = antecedents.end();
      for (const_iterator aiter = antecedents.begin(); aiter != aiter_end; ++ aiter)
	aiter->frontier(iter);
    }

  private:
    size_type __depth(const TreeRule& tree, const size_type depth) const
    {
      size_type max_depth = depth;
      
      const_iterator aiter_end = tree.end();
      for (const_iterator aiter = tree.begin(); aiter != aiter_end; ++ aiter)
	max_depth = utils::bithack::max(max_depth, __depth(*aiter, depth + 1));
      
      return max_depth;
    }

    template <typename Path>
    void __hyperpath(const TreeRule& tree, Path& path, const size_type max_depth, const size_type depth) const
    {
      path[depth].push_back(tree.label);
      
      if (tree.empty()) {
	for (size_type i = depth; i != max_depth; ++ i) {
	  path[i + 1].push_back(vocab_type::EPSILON);
	  path[i + 1].push_back(vocab_type::NONE);
	}
      } else {
	const_iterator aiter_end = tree.end();
	for (const_iterator aiter = tree.begin(); aiter != aiter_end; ++ aiter)
	  __hyperpath(*aiter, path, max_depth, depth + 1);
	path[depth + 1].push_back(vocab_type::NONE);
      }
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
  size_t __hash_value_tree_rule(TreeRule const& x, size_t seed)
  {
    for (TreeRule::const_iterator aiter = x.begin(); aiter != x.end(); ++ aiter)
      seed = __hash_value_tree_rule(*aiter, seed);
    
    return utils::hashmurmur<size_t>()(x.label.id(), seed);
  }
  
  inline
  size_t hash_value(TreeRule const& x)
  {
    return __hash_value_tree_rule(x, 0);
  }

  void sort(TreeRule& x, TreeRule& y);
  
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
