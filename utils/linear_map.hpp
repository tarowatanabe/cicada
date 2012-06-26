// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__LINEAR_MAP__HPP__
#define __UTILS__LINEAR_MAP__HPP__ 1

#include <utility>
#include <functional>
#include <algorithm>
#include <vector>

#include <utils/simple_vector.hpp>

#include <boost/functional/hash.hpp>

namespace utils
{
  
  template <typename Key, typename Data, typename Predicate>
  struct __linear_map_predicate_impl
  {
    struct predicate : public Predicate
    {
      bool operator()(const std::pair<Key, Data>& x, const std::pair<Key, Data>& y) const
      {
	return Predicate::operator()(x.first, y.first);
      }
      bool operator()(const std::pair<Key, Data>& x, const Key& y) const
      {
	return Predicate::operator()(x.first, y);
      }
      bool operator()(const Key& x, const std::pair<Key, Data>& y) const
      {
	return Predicate::operator()(x, y.first);
      }
      bool operator()(const Key& x, const Key& y) const
      {
	return Predicate::operator()(x, y);
      }
    };
  };

  template <typename Key,
	    typename Data,
	    typename Hash=boost::hash<Key>,
	    typename Predicate=std::equal_to<Key>,
	    typename Alloc=std::allocator<std::pair<Key, Data> > >
  class linear_map : public __linear_map_predicate_impl<Key,Data,Predicate>::predicate, public Hash
  {
  public:
    typedef Key                              key_type;
    typedef Data                             mapped_type;
    typedef Data                             data_type;
    typedef std::pair<key_type, mapped_type> value_type;
    
    typedef Hash      hash_type;
    
  private:
    typedef typename __linear_map_predicate_impl<Key,Data,Predicate>::predicate predicate_type;
    
    typedef typename Alloc::template rebind<value_type>::other alloc_type;
    typedef utils::simple_vector<value_type, alloc_type> container_type;

  public:
    typedef typename container_type::pointer pointer;
    typedef typename container_type::const_pointer const_pointer;
    typedef typename container_type::reference reference;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;
    typedef typename container_type::reverse_iterator reverse_iterator;
    typedef typename container_type::const_reverse_iterator const_reverse_iterator;
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;

  public:
    linear_map() : container() {}
    template <typename _InputIterator>
    linear_map(_InputIterator first, _InputIterator last) : container()
    {
      for (/**/; first != last; ++ first)
	insert(*first);
    }
    
    void swap(linear_map& x)
    {
      container.swap(x.container);
      std::swap(static_cast<hash_type&>(*this), static_cast<hash_type&>(x));
      std::swap(static_cast<predicate_type&>(*this), static_cast<predicate_type&>(x));
    }


    iterator begin() { return container.begin(); }
    const_iterator begin() const { return container.begin(); }
    iterator end() { return container.end(); }
    const_iterator end() const { return container.end(); }
    reverse_iterator rbegin() { return container.rbegin(); }
    const_reverse_iterator rbegin() const { return container.rbegin(); }
    reverse_iterator rend() { return container.rend(); }
    const_reverse_iterator rend() const { return container.rend(); }
    
    // non-standard!
    reference front() { return container.front(); }
    const_reference front() const { return container.front(); }
    reference back() { return container.back(); }
    const_reference back() const { return container.back(); }
    
    bool empty() const { return container.empty(); }
    size_type size() const { return container.size(); }
    
    void erase(iterator position)
    {
      container.erase(position);
    }
    size_type erase(const key_type& x)
    {
      iterator iter_end = end();
      for (iterator iter = begin(); iter != iter_end; ++ iter)
	if (predicate_type::operator()(*iter, x)) {
	  container.erase(iter);
	  return 1;
	}
      return 0;
    }
    void erase(iterator first, iterator last)
    {
      container.erase(first, last);
    }
    void clear()
    {
      container.clear();
    }
    
    Data& operator[](const key_type& k)
    {
      iterator iter_end = end();
      for (iterator iter = begin(); iter != iter_end; ++ iter)
	if (predicate_type::operator()(*iter, k))
	  return iter->second;
      
      return container.insert(iter_end, std::make_pair(k, data_type()))->second;
    }
    
    std::pair<iterator, bool> insert(const value_type& x)
    {
      iterator iter_end = end();
      for (iterator iter = begin(); iter != iter_end; ++ iter)
	if (predicate_type::operator()(*iter, x))
	  return std::make_pair(iter, false);
      
      return std::make_pair(container.insert(iter_end, x), true);
    }
    
    iterator insert(iterator position, const value_type& x)
    {
      return insert(x).first;
    }
    template <typename InputIterator>
    void insert(InputIterator first, InputIterator last)
    {
      for (/**/; first != last; ++ first)
	insert(*first);
    }
    
    iterator find(const key_type& x)
    {
      iterator iter_end = end();
      for (iterator iter = begin(); iter != iter_end; ++ iter)
	if (predicate_type::operator()(*iter, x))
	  return iter;
      return iter_end;
    }
    const_iterator find(const key_type& x) const
    {
      const_iterator iter_end = end();
      for (const_iterator iter = begin(); iter != iter_end; ++ iter)
	if (predicate_type::operator()(*iter, x))
	  return iter;
      return iter_end;
    }
    
  public:
    
    template <typename K, typename D, typename H, typename P, typename A>
    friend
    bool operator==(const linear_map<K,D,H,P,A>& x, const linear_map<K,D,H,P,A>& y);
    
    template <typename K, typename D, typename H, typename P, typename A>
    friend
    bool operator!=(const linear_map<K,D,H,P,A>& x, const linear_map<K,D,H,P,A>& y);
    
  private:
    container_type container;
  };

  template <typename K, typename D, typename H, typename P, typename A>
  inline
  bool operator==(const linear_map<K,D,H,P,A>& x, const linear_map<K,D,H,P,A>& y)
  {
    if (x.size() != y.size())
      return false;
    else {
      typedef linear_map<K,D,H,P,A> container_type;
      
      typename container_type::const_iterator iter_end = y.end();
      for (typename container_type::const_iterator iter = y.begin(); iter != iter_end; ++ iter)
	if (x.find(iter->first) == x.end())
	  return false;
    }
    return true;
  }
  
  template <typename K, typename D, typename H, typename P, typename A>
  inline
  bool operator!=(const linear_map<K,D,H,P,A>& x, const linear_map<K,D,H,P,A>& y)
  {
    return ! (x == y);
  }
};

namespace std
{
  template <typename K, typename D, typename H, typename P, typename A>
  inline
  void swap(utils::linear_map<K,D,H,P,A>& x, 
	    utils::linear_map<K,D,H,P,A>& y)
  {
    x.swap(y);
  }
};

#endif
