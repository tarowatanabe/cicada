// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__VECTOR2__HPP__
#define __UTILS__VECTOR2__HPP__ 1

#include <vector>
#include <boost/swap.hpp>

namespace utils
{
  
  template <typename Iterator, typename Reference>
  struct __vector2_view
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Iterator  iterator;
    typedef Reference reference;

    typedef __vector2_view self_type;
    
    __vector2_view(iterator __iter, size_type __size) : iter(__iter), size(__size) {}
    
    reference operator*() const { return *iter; }
    iterator  operator->() const { return iter; }
    reference operator[](difference_type __n) const { return *(iter + __n); }
    
    self_type& operator++() { iter += size; return *this; }
    self_type& operator--() { iter -= size; return *this; }
    self_type  operator++(int) { self_type __tmp = *this; ++ *this; return __tmp; }
    self_type  operator--(int) { self_type __tmp = *this; -- *this; return __tmp; }
    self_type& operator+=(difference_type __n) { iter += __n * size; return *this; }
    self_type& operator-=(difference_type __n) { iter -= __n * size; return *this; }
    self_type operator+(difference_type __n) const { self_type __tmp = *this; return __tmp += __n; }
    self_type operator-(difference_type __n) const { self_type __tmp = *this; return __tmp -= __n; }
    
    
    iterator iter;
    size_type size;
  };

  template <typename _Tp, typename _Alloc=std::allocator<_Tp> >
  class vector2
  {
  private:
    typedef std::vector<_Tp, _Alloc> base_type;
    
  public:
    typedef typename base_type::value_type      value_type;
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::difference_type difference_type;
    
    typedef typename base_type::const_iterator  const_iterator;
    typedef typename base_type::iterator        iterator;
    typedef typename base_type::const_reference const_reference;
    typedef typename base_type::reference       reference;
    
    typedef __vector2_view<const_iterator, const_reference> const_view_type;
    typedef __vector2_view<iterator, reference>             view_type;
    
  public:
    vector2()
      : __base(), __size1(), __size2() {}
    vector2(size_type __s1, size_type __s2)
      : __base(__s1 * __s2), __size1(__s1), __size2(__s2) {}
    vector2(size_type __s1, size_type __s2, const value_type& value)
      : __base(__s1 * __s2, value), __size1(__s1), __size2(__s2) {}
    
    bool empty() const { return __base.empty(); }
    size_type size1() const { return __size1; }
    size_type size2() const { return __size2; }
    
    void reserve(size_type __s1, size_type __s2)
    {
      __base.reserve(__s1 * __s2);
    }
    
    void resize(size_type __s1, size_type __s2)
    {
      __size1 = __s1;
      __size2 = __s2;
      __base.resize(__s1 * __s2);
    }
    void resize(size_type __s1, size_type __s2, const value_type& value)
    {
      __size1 = __s1;
      __size2 = __s2;
      __base.resize(__s1 * __s2, value);
    }
    void clear()
    {
      __base.clear();
      __size1 = 0;
      __size2 = 0;
    }
    
    void swap(vector2& x)
    {
      __base.swap(x.__base);
      boost::swap(__size1, x.__size1);
      boost::swap(__size2, x.__size2);
    }

    void swap(size_type pos1, size_type pos2)
    {
      if (pos1 == pos2) return;
      
      for (size_type j = 0; j < __size2; ++ j)
	boost::swap(operator()(pos1, j), operator()(pos2, j));
      
      for (size_type i = 0; i < __size1; ++ i)
	boost::swap(operator()(i, pos1), operator()(i, pos2));
    }
    
    inline const_view_type operator[](size_type pos1) const { return const_view_type(begin(pos1), __size2); }
    inline       view_type operator[](size_type pos1)       { return view_type(begin(pos1), __size2); }
    
    inline const_reference operator()(size_type pos1, size_type pos2) const { return __base[pos1 * __size2 + pos2]; }
    inline       reference operator()(size_type pos1, size_type pos2)       { return __base[pos1 * __size2 + pos2]; }
    
    inline const_iterator begin() const { return __base.begin(); }
    inline       iterator begin()       { return __base.begin(); }
    
    inline const_iterator end() const { return __base.end(); }
    inline       iterator end()       { return __base.end(); }
    
    inline const_iterator begin(size_type pos) const { return __base.begin() + pos * __size2; }
    inline       iterator begin(size_type pos)       { return __base.begin() + pos * __size2; }
    
    inline const_iterator end(size_type pos) const { return __base.begin() + (pos + 1) * __size2; }
    inline       iterator end(size_type pos)       { return __base.begin() + (pos + 1) * __size2; }
    
  private:
    base_type __base;
    size_type __size1;
    size_type __size2;
  };
  
};

namespace std
{
  template <typename Tp, typename Alloc>
  inline
  void swap(utils::vector2<Tp, Alloc>& x, utils::vector2<Tp, Alloc>& y)
  {
    x.swap(x);
  }
};
#endif
