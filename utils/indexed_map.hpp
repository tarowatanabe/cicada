// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__INDEXED_MAP__HPP__
#define __UTILS__INDEXED_MAP__HPP__ 1

#include <utils/indexed_hashtable.h>

#include <boost/functional/hash/hash.hpp>

namespace utils
{
  
  template <typename Key,
	    typename Data,
	    typename Hash=boost::hash<Key>,
	    typename Equal=std::equal_to<Key>,
	    typename Alloc=std::allocator<std::pair<const Key, Data> > >
  class indexed_map 
  {
  public:
    typedef Key                                  key_type;
    typedef Data                                 mapped_type;
    typedef std::pair<const key_type, mapped_type> value_type;
    
  private:
    struct extract_key
    {
      const Key& operator()(const value_type& x) const { return x.first; }
      const Key& operator()(value_type& x) const { return x.first; }
    };
    typedef indexed_hashtable<key_type, value_type, extract_key, Hash, Equal, Alloc> impl_type;

  public:
    typedef typename impl_type::size_type  size_type;
    typedef typename impl_type::index_type index_type;
    
    typedef typename impl_type::iterator       iterator;
    typedef typename impl_type::const_iterator const_iterator;
    
    typedef mapped_type&      reference;
    typedef const mapped_type& const_reference;

  public:
    indexed_map(const size_type __size=8, const Hash& __hash=Hash(), const Equal& __equal=Equal())
      : impl(__size, __hash, __equal) {}

  public:
    void assign(const indexed_map& x) { impl.assign(x.impl); }
    void swap(indexed_map& x) { impl.swap(x.impl); }
    
  public:
    
    inline const value_type& operator[](index_type x) const { return impl[x]; }
    inline       value_type& operator[](index_type x)       { return impl[x]; }
    
    const_iterator begin() const { return impl.begin(); }
    const_iterator end() const { return impl.end(); }
    
    bool empty() const { return impl.empty(); }
    size_type size() const { return impl.size(); }
    void clear() { impl.clear(); }
    
    std::pair<iterator, bool> insert(const value_type& x) { return impl.insert(x); }
    const_iterator find(const key_type& x) const { return impl.find(x); }
    
  private:
    impl_type impl;
  };
  
};

namespace std
{
  template <typename Key, typename Data, typename Hash, typename Equal, typename Alloc>
  inline
  void swap(utils::indexed_map<Key,Data, Hash,Equal,Alloc>& x,
	    utils::indexed_map<Key,Data, Hash,Equal,Alloc>& y)
  {
    x.swap(y);
  }
  
};

#endif



