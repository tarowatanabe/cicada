// -*- mode: c++ -*-

#ifndef __UTILS__DENSE_SET__H__
#define __UTILS__DENSE_SET__H__ 1

#include <utility>
#include <algorithm>
#include <functional>

#include <utils/dense_hashtable.hpp>

#include <boost/functional/hash.hpp>

namespace utils
{
  
  template <typename Key,
	    typename Hash=boost::hash<Key>,
	    typename Equal=std::equal_to<Key>,
	    typename Alloc=std::allocator<std::pair<const Key, Data> > >
  class dense_set
  {
  public:
    typedef Key  key_type;
    typedef Key  data_type;
    typedef Key  value_type;
    
  private:
    struct extract_key
    {
      const key_type& operator()(const value_type& x) const
      {
	return x;
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
    
    inline const_iterator find(const key_type& x) const { return hashtable.find(x); }
    inline       iterator find(const key_type& x)       { return hashtable.find(x); }
    
    bool empty() const { return hashtable.empty(); }
    size_type size() const { return hashtable.size(); }
    
    
    void clear() { hashtable.clear(); }
    void swap(dense_set& x) { hashtable.swap(x.hashtable); }
    
    void insert(const value_type& x) { return hashtable.erase(x); }
    
    void erase(const key_type& x) { hashtable.erase(x); }
    void erase(iterator x) { hashtable.erase(x); }
    
  private:
    hashtable_type hashtable;
  };
  
};

namespace std
{
  template <typename K, typename H, typename E, typename A>
  inline
  void swap(utils::dense_set<K,H,E,A>& x,
	    utils::dense_set<K,H,E,A>& y)
  {
    x.swap(y);
  }
};



#endif
