// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__LINEAR_SET__HPP__
#define __UTILS__LINEAR_SET__HPP__ 1

#include <utility>
#include <functional>
#include <algorithm>
#include <vector>

#include <utils/simple_vector.hpp>

#include <boost/functional/hash.hpp>

namespace utils
{
  template <typename Tp,
	    typename Hash=boost::hash<Tp>,
	    typename Predicate=std::equal_to<Tp>,
	    typename Alloc=std::allocator<Tp > >
  class linear_set : public Hash,
		     public Predicate
  {
  public:
    typedef Tp        value_type;
    typedef Hash      hash_type;
    typedef Predicate predicate_type;
    
  private:
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
    linear_set() : container() {}
    template <typename _InputIterator>
    linear_set(_InputIterator first, _InputIterator last) : container()
    {
      for (/**/; first != last; ++ first)
	insert(*first);
    }
    
    void swap(linear_set& x)
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
    size_type erase(const value_type& x)
    {
      iterator iter_end = end();
      for (iterator iter = begin(); iter != iter_end; ++ iter)
	if (predicate()(*iter, x)) {
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
    
    std::pair<iterator, bool> insert(const value_type& x)
    {
      iterator iter_end = end();
      for (iterator iter = begin(); iter != iter_end; ++ iter)
	if (predicate()(*iter, x))
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
    
    iterator find(const value_type& x)
    {
      iterator iter_end = end();
      for (iterator iter = begin(); iter != iter_end; ++ iter)
	if (predicate()(*iter, x))
	  return iter;
      return iter_end;
    }
    
    const_iterator find(const value_type& x) const
    {
      const_iterator iter_end = end();
      for (const_iterator iter = begin(); iter != iter_end; ++ iter)
	if (predicate()(*iter, x))
	  return iter;
      return iter_end;
    }

  public:
    const hash_type& hasher() const { return static_cast<const hash_type&>(*this); }
    hash_type& hasher() { return static_cast<hash_type&>(*this); }

    const predicate_type& predicate() const { return static_cast<const predicate_type&>(*this); }
    predicate_type& predicate() { return static_cast<predicate_type&>(*this); }

  public:
    template <typename T, typename H, typename P, typename A>
    friend
    bool operator==(const linear_set<T,H,P,A>& x, const linear_set<T,H,P,A>& y);
    
    template <typename T, typename H, typename P, typename A>
    friend
    bool operator!=(const linear_set<T,H,P,A>& x, const linear_set<T,H,P,A>& y);
    
  private:
    container_type container;
  };

  template <typename T, typename H, typename P, typename A>
  inline
  bool operator==(const linear_set<T,H,P,A>& x, const linear_set<T,H,P,A>& y)
  {
    if (x.size() != y.size())
      return false;
    else if (std::equal(x.begin(), x.end(), y.begin(), x.predicate()))
      return true;
    else {
      typedef linear_set<T,H,P,A> container_type;
      
      typename container_type::const_iterator iter_end = y.end();
      for (typename container_type::const_iterator iter = y.begin(); iter != iter_end; ++ iter)
	if (x.find(*iter) == x.end())
	  return false;
    }
    return true;
  }
  
  template <typename T, typename H, typename P, typename A>
  inline
  bool operator!=(const linear_set<T,H,P,A>& x, const linear_set<T,H,P,A>& y)
  {
    return ! (x == y);
  }

};

namespace std
{
  template <typename T, typename H, typename P, typename A>
  inline
  void swap(utils::linear_set<T,H,P,A>& x, 
	    utils::linear_set<T,H,P,A>& y)
  {
    x.swap(y);
  }
};

#endif
