// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__STD_HEAP__HPP__
#define __UTILS__STD_HEAP__HPP__ 1

#include <vector>
#include <queue>
#include <functional>

#include <utils/bithack.hpp>

namespace utils
{
  // priority_queue like interface...
  template <typename _Tp, typename _Sequence=std::vector<_Tp>, typename _Compare=std::less<typename _Sequence::value_type> >
  class std_heap : public _Compare
  {
  public:
    typedef typename _Sequence::value_type                value_type;
    typedef typename _Sequence::reference                 reference;
    typedef typename _Sequence::const_reference           const_reference;
    typedef typename _Sequence::size_type                 size_type;
    typedef          _Sequence                            container_type;

    typedef typename _Sequence::const_iterator const_iterator;
    typedef typename _Sequence::iterator       iterator;

  protected:
    _Sequence c;
    
    
  public:
    std_heap(size_type x) : c(x) { clear(); }
    std_heap() {  }
    std_heap(size_type x, const _Compare& __comp) : _Compare(__comp), c(x) { clear(); }
    std_heap(const _Compare& __comp) : _Compare(__comp) { }
    
    bool empty() const { return c.empty(); }
    size_type size() const { return c.size(); }
    
    void clear() { c.clear(); }
    
    void reserve(size_type x) { c.reserve(x); }
    
    // a conventional interface...
    const_reference top()
    {
      return c.front();
    }

    void push(const value_type& x)
    {
      c.push_back(x);
      std::push_heap(c.begin(), c.end(), static_cast<_Compare&>(*this));
    }
    
    void pop()
    {
      std::pop_heap(c.begin(), c.end(), static_cast<_Compare&>(*this));
      c.pop_back();
    }

    void push_back(const value_tyep& x)
    {
      c.push_back(x);
    }
    
    void make_heap()
    {
      std::make_heap(c.begin(), c.end(), static_cast<_Compare&>(*this));
    }
    
    
    const_iterator begin() const { return c.begin(); }
    const_iterator end() const { return c.end(); }
    
    void swap(std_heap& x)
    {
      std::swap(static_cast<_Compare&>(*this), static_cast<_Compare&>(x));
      c.swap(x.c);
    }

  };
  
};

namespace std
{
  template <typename T, typename S, typename C>
  inline
  void swap(utils::std_heap<T,S,C>& x,
	    utils::std_heap<T,S,C>& y)
  {
    x.swap(y);
  }
};


#endif
