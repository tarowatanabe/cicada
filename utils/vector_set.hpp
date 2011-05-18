// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__VECTOR_SET__HPP__
#define __UTILS__VECTOR_SET__HPP__ 1

#include <utility>
#include <functional>
#include <algorithm>
#include <vector>
#include <utils/simple_vector.hpp>

namespace utils
{
  
  
  template <typename Tp, typename Compare=std::less<Tp>, typename Alloc=std::allocator<Tp > >
  class vector_set : public Compare
  {
  public:
    typedef Tp      value_type;
    typedef Compare compare_type;
    
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
    vector_set() : container() {}
    template <typename _InputIterator>
    vector_set(_InputIterator first, _InputIterator last) : container(first, last)
    {
      stable_sort();
      unique();
    }
    
    void swap(vector_set& x)
    {
      container.swap(x.container);
      std::swap(static_cast<compare_type&>(*this), static_cast<compare_type&>(x));
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
      std::pair<iterator, iterator> range = equal_range(x);
      if (range.first != range.second)
	erase(range.first, range.second);
      return range.second - range.first;
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
      iterator iter = lower_bound(x);
      if (iter != container.end() && (! compare_type::operator()(x, *iter)))
	return std::make_pair(iter, false);
      return std::make_pair(container.insert(iter, x), true);
    }
    
    iterator insert(iterator position, const value_type& x)
    {
      if (position == end()) {
	if (! empty() && compare_type::operator()(*(position - 1), x))
	  return container.insert(position, x);
	else
	  return insert(x).first;
      } else if (position == begin()) {
	if (! empty() && compare_type::operator()(x, *position))
	  return container.insert(position, x);
	else
	  return insert(x).first;
      } else {
	if (compare_type::operator()(*(position - 1), x) && compare_type::operator()(x, *position))
	  return container.insert(position, x);
	else
	  return insert(x).first;
      }
    }
    template <typename InputIterator>
    void insert(InputIterator first, InputIterator last)
    {
      container.insert(container.end(), first, last);
      stable_sort();
      unique();
    }
    
    iterator find(const value_type& x)
    {
      iterator iter = lower_bound(x);
      if (iter != container.end() && (! compare_type::operator()(x, *iter)))
	return iter;
      return end();
    }
    const_iterator find(const value_type& x) const
    {
      const_iterator iter = lower_bound(x);
      if (iter != container.end() && (! compare_type::operator()(x, *iter)))
	return iter;
      return end();
    }

    iterator lower_bound(const value_type& x)
    {
      if (container.size() < 16) {
	iterator iter = container.begin();
	iterator end = container.end();
	for (/**/; iter != end && static_cast<const compare_type&>(*this)(*iter, x); ++ iter);
	return iter;
      } else
	return std::lower_bound(container.begin(), container.end(), x, static_cast<const compare_type&>(*this));
    }
    const_iterator lower_bound(const value_type& x) const
    {
      if (container.size() < 16) {
	const_iterator iter = container.begin();
	const_iterator end = container.end();
	for (/**/; iter != end && static_cast<const compare_type&>(*this)(*iter, x); ++ iter);
	return iter;
      } else
	return std::lower_bound(container.begin(), container.end(), x, static_cast<const compare_type&>(*this));
    }
  
    iterator upper_bound(const value_type& x)
    {
      if (container.size() < 16) {
	iterator iter = container.begin();
	iterator end = container.end();
	for (/**/; iter != end && static_cast<const compare_type&>(*this)(*iter, x); ++ iter);
	for (/**/; iter != end && ! static_cast<const compare_type&>(*this)(x, *iter); ++ iter);
	return iter;
      } else
	return std::upper_bound(container.begin(), container.end(), x, static_cast<const compare_type&>(*this));  
    }
    const_iterator upper_bound(const value_type& x) const
    {
      if (container.size() < 16) {
	const_iterator iter = container.begin();
	const_iterator end = container.end();
	for (/**/; iter != end && static_cast<const compare_type&>(*this)(*iter, x); ++ iter);
	for (/**/; iter != end && ! static_cast<const compare_type&>(*this)(x, *iter); ++ iter);
	return iter;
      } else
	return std::upper_bound(container.begin(), container.end(), x, static_cast<const compare_type&>(*this));
    }
  
    std::pair<iterator, iterator> equal_range(const value_type& x)
    {
      if (container.size() < 16) {
	iterator iter = container.begin();
	iterator end = container.end();
	for (/**/; iter != end && static_cast<const compare_type&>(*this)(*iter, x); ++ iter);
	iterator first = iter;
	for (/**/; iter != end && ! static_cast<const compare_type&>(*this)(x, *iter); ++ iter);
	return std::make_pair(first, iter);
      } else
	return std::equal_range(container.begin(), container.end(), x, static_cast<const compare_type&>(*this));
    }
    std::pair<const_iterator, const_iterator> equal_range(const value_type& x) const
    {
      if (container.size() < 16) {
	const_iterator iter = container.begin();
	const_iterator end = container.end();
	for (/**/; iter != end && static_cast<const compare_type&>(*this)(*iter, x); ++ iter);
	const_iterator first = iter;
	for (/**/; iter != end && ! static_cast<const compare_type&>(*this)(x, *iter); ++ iter);
	return std::make_pair(first, iter);
      } else
	return std::equal_range(container.begin(), container.end(), x, static_cast<const compare_type&>(*this));
    }

  public:
    
    template <typename T, typename C, typename A>
    friend
    bool operator==(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y);
    
    template <typename T, typename C, typename A>
    friend
    bool operator!=(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y);
    
    template <typename T, typename C, typename A>
    friend
    bool operator<(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y);
    
    template <typename T, typename C, typename A>
    friend
    bool operator>(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y);
    
    template <typename T, typename C, typename A>
    friend
    bool operator<=(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y);
    
    template <typename T, typename C, typename A>
    friend
    bool operator>=(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y);
    
  private:
    void stable_sort()
    {
      if (container.size() > static_cast<size_type>(1))
	std::stable_sort(container.begin(), container.end(), static_cast<const compare_type&>(*this));
    }
    
    void unique()
    {
      // we assume that we are already sorted
      
      if (container.size() <= static_cast<size_type>(1)) return;
      
      iterator iter1 = container.begin();
      iterator iter2 = container.begin();
      while (iter2 != container.end()) {
	iterator bound = std::upper_bound(iter2, container.end(), *iter2, static_cast<const compare_type&>(*this));
	if (iter1 != iter2)
	  *iter1 = *iter2;
	++ iter1;
	iter2 = bound;
      }
      
      container_type container_new(container.begin(), iter1);
      container.swap(container_new);
    }
    
  private:
    container_type container;
  };

  template <typename T, typename C, typename A>
  inline
  bool operator==(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y)
  {
    return x.container == y.container;
  }
  
  template <typename T, typename C, typename A>
  inline
  bool operator!=(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y)
  {
    return x.container != y.container;
  }

  template <typename T, typename C, typename A>
  inline
  bool operator<(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y)
  {
    return x.container < y.container;
  }
  
  template <typename T, typename C, typename A>
  inline
  bool operator>(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y)
  {
    return x.container > y.container;
  }

  template <typename T, typename C, typename A>
  inline
  bool operator<=(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y)
  {
    return x.container <= y.container;
  }
  
  template <typename T, typename C, typename A>
  inline
  bool operator>=(const vector_set<T,C,A>& x, const vector_set<T,C,A>& y)
  {
    return x.container >= y.container;
  }

};

namespace std
{
  template <typename T, typename C, typename A>
  inline
  void swap(utils::vector_set<T,C,A>& x, 
	    utils::vector_set<T,C,A>& y)
  {
    x.swap(y);
  }
};

#endif
