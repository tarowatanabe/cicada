// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__DENSE_MAP__H__
#define __UTILS__DENSE_MAP__H__ 1

#include <utility>
#include <algorithm>
#include <functional>

#include <utils/dense_hashtable.hpp>

#include <boost/functional/hash.hpp>

namespace utils
{
  
  template <typename Key,
	    typename Data,
	    typename Hash=boost::hash<Key>,
	    typename Equal=std::equal_to<Key>,
	    typename Alloc=std::allocator<std::pair<const Key, Data> > >
  class dense_map
  {
  public:
    typedef Key  key_type;
    typedef Data data_type;
    typedef Data mapped_type;
    typedef std::pair<const key_type, data_type> value_type;
    
  private:
    struct extract_key
    {
      const key_type& operator()(const value_type& x) const
      {
	return x.first;
      }
    };
    typedef utils::dense_hashtable<key_type, value_type, extract_key, Hash, Equal, Alloc> hashtable_type;

  public:
    typedef typename hashtable_type::size_type       size_type;
    typedef typename hashtable_type::difference_type difference_type;
    
    typedef typename hashtable_type::const_iterator const_iterator;
    typedef typename hashtable_type::iterator       iterator;
    
    typedef       value_type& reference;
    typedef const value_type& const_reference;

  public:
    
    inline const_iterator begin() const { return hashtable.begin(); }
    inline       iterator begin()       { return hashtable.begin(); }
    
    inline const_iterator end() const { return hashtable.end(); }
    inline       iterator end()       { return hashtable.end(); }
    
    inline const_iterator find(const key_type& x) const { return hashtable.find(value_type(x, data_type())); }
    inline       iterator find(const key_type& x)       { return hashtable.find(value_type(x, data_type())); }
    
    bool empty() const { return hashtable.empty(); }
    size_type size() const { return hashtable.size(); }
    
    data_type& operator[](const key_type& x)
    {
      return hashtable.insert(std::make_pair(x, data_type())).first->second;
    }
    
    void clear() { hashtable.clear(); }
    void swap(dense_map& x) { hashtable.swap(x.hashtable); }
    
    void insert(const value_type& x) { return hashtable.erase(x); }
    
    void erase(const key_type& x) { hashtable.erase(value_type(x, data_type())); }
    void erase(iterator x) { hashtable.erase(x); }
    
  private:
    hashtable_type hashtable;
  };
  
};

namespace std
{
  template <typename K, typename D, typename H, typename E, typename A>
  inline
  void swap(utils::dense_map<K,D,H,E,A>& x,
	    utils::dense_map<K,D,H,E,A>& y)
  {
    x.swap(y);
  }
};



#endif
