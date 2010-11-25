// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__INDEXED_SET__HPP__
#define __UTILS__INDEXED_SET__HPP__ 1

#include <utils/indexed_hashtable.hpp>

#include <boost/functional/hash/hash.hpp>

namespace utils
{
  template <typename Tp, typename Hash=boost::hash<Tp>, typename Equal=std::equal_to<Tp>, typename Alloc=std::allocator<Tp> >
  class indexed_set
  {
    private:
    struct extract_key
    {
      const Tp& operator()(const Tp& x) const { return x; }
      Tp& operator()(Tp& x) const { return x; }
    };
    typedef indexed_hashtable<Tp, Tp, extract_key, Hash, Equal, Alloc> impl_type;
    
  public:
    typedef Tp value_type;
    
    typedef typename impl_type::size_type  size_type;
    typedef typename impl_type::index_type index_type;
    
    typedef typename impl_type::const_iterator iterator;
    typedef typename impl_type::const_iterator const_iterator;
    typedef typename impl_type::const_reference reference;
    typedef typename impl_type::const_reference const_reference;

  public:
    indexed_set(const size_type __size=8, const Hash& __hash=Hash(), const Equal& __equal=Equal())
      : impl(__size, __hash, __equal) {}
    
    
  public:
    void assign(const indexed_set& x) { impl.assign(x.impl); }
    void swap(indexed_set& x) { impl.swap(x.impl); }
    
  public:
    
    const Tp& operator[](index_type x) const { return impl[x]; }

    const_iterator begin() const { return impl.begin(); }
    const_iterator end() const { return impl.end(); }
    
  public:
    bool empty() const { return impl.empty(); }
    size_type size() const { return impl.size(); }
    void clear() { impl.clear(); }
    
    std::pair<iterator, bool> insert(const value_type& x) { return impl.insert(x); }
    const_iterator find(const value_type& x) const { return impl.find(x); }
    
  private:
    impl_type impl;    
  };
  
};

namespace std
{
  template <typename Tp, typename Hash, typename Equal, typename Alloc>
  inline
  void swap(utils::indexed_set<Tp,Hash,Equal,Alloc>& x,
	    utils::indexed_set<Tp,Hash,Equal,Alloc>& y)
  {
    x.swap(y);
  }
  
};

#endif
