// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__VECTOR2_ALIGNED__HPP__
#define __UTILS__VECTOR2_ALIGNED__HPP__ 1

#include <vector>

namespace utils
{
  
  template <typename Iterator, typename Reference>
  struct __vector2_aligned_view
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Iterator  iterator;
    typedef Reference reference;

    typedef __vector2_aligned_view self_type;
    
    __vector2_aligned_view(iterator __iter, size_type __size) : iter(__iter), size(__size) {}
    
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
  class vector2_aligned
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
    
    typedef __vector2_aligned_view<const_iterator, const_reference> const_view_type;
    typedef __vector2_aligned_view<iterator, reference>             view_type;
    
  public:
    vector2_aligned()
      : __base(), __size1(0), __size2(0), __size1_aligned(0), __size2_aligned(0) {}
    vector2_aligned(size_type __s1, size_type __s2)
      : __base(__vector_size(__s1, __s2)),
	__size1(__s1), __size2(__s2),
	__size1_aligned(__s1), __size2_aligned(__aligned_size(__s2)) {}
    vector2_aligned(size_type __s1, size_type __s2, const value_type& value)
      : __base(__vector_size(__s1, __s2), value),
	__size1(__s1), __size2(__s2),
	__size1_aligned(__s1), __size2_aligned(__aligned_size(__s2)) {}
    
    bool empty() const { return __base.empty(); }
    size_type size1() const { return __size1; }
    size_type size2() const { return __size2; }
    
    void reserve(size_type __s1, size_type __s2)
    {
      __base.reserve(__vector_size(__s1, __s2));
    }
    
    void resize(size_type __s1, size_type __s2)
    {
      __size1 = __s1;
      __size2 = __s2;
      
      __size1_aligned = __size1;
      __size2_aligned = __aligned_size(__size2);
      
      __base.resize(__size1_aligned * __size2_aligned);
    }
    void resize(size_type __s1, size_type __s2, const value_type& value)
    {
      __size1 = __s1;
      __size2 = __s2;
      
      __size1_aligned = __size1;
      __size2_aligned = __aligned_size(__size2);
      
      __base.resize(__size1_aligned * __size2_aligned, value);
    }
    void clear()
    {
      __base.clear();
      __size1 = 0;
      __size2 = 0;
      __size1_aligned = 0;
      __size2_aligned = 0;
    }
    void swap(vector2_aligned& x)
    {
      __base.swap(x.__base);
      
      std::swap(__size1, x.__size1);
      std::swap(__size2, x.__size2);
      std::swap(__size1_aligned, x.__size1_aligned);
      std::swap(__size2_aligned, x.__size2_aligned);
    }

    void swap(size_type pos1, size_type pos2)
    {
      using namespace std;
      
      if (pos1 == pos2) return;
      
      for (size_type j = 0; j < __size2; ++ j)
        swap(operator()(pos1, j), operator()(pos2, j));
      
      for (size_type i = 0; i < __size1; ++ i)
        swap(operator()(i, pos1), operator()(i, pos2));
    }
    
    inline const_view_type operator[](size_type pos1) const { return const_view_type(begin(pos1), __size2_aligned); }
    inline       view_type operator[](size_type pos1)       { return view_type(begin(pos1), __size2_aligned); }
    
    inline const_reference operator()(size_type pos1, size_type pos2) const { return __base[pos1 * __size2_aligned + pos2]; }
    inline       reference operator()(size_type pos1, size_type pos2)       { return __base[pos1 * __size2_aligned + pos2]; }
    
    inline const_iterator begin() const { return __base.begin(); }
    inline       iterator begin()       { return __base.begin(); }
    
    inline const_iterator end() const { return __base.end(); }
    inline       iterator end()       { return __base.end(); }
    
    inline const_iterator begin(size_type pos) const { return __base.begin() + pos * __size2_aligned; }
    inline       iterator begin(size_type pos)       { return __base.begin() + pos * __size2_aligned; }
    
    inline const_iterator end(size_type pos) const { return __base.begin() + (pos + 1) * __size2_aligned; }
    inline       iterator end(size_type pos)       { return __base.begin() + (pos + 1) * __size2_aligned; }
    
  private:
    size_type __vector_size(size_type __s1, size_type __s2) const
    {
      return __s1 * __aligned_size(__s2);
    }
    
    size_type __aligned_size(size_type size) const
    {
      static const int alignment_size = 16 / sizeof(value_type);

      return (size + alignment_size - 1) & (size_type(-alignment_size));
    }
    
  private:
    base_type __base;
    
    size_type __size1;
    size_type __size2;
    
    size_type __size1_aligned; // this will be the same as the size1
    size_type __size2_aligned; // this will be the 16-byte alignment adjusted size2
  };
};

namespace std
{
  template <typename Tp, typename Alloc>
  inline
  void swap(utils::vector2_aligned<Tp, Alloc>& x, utils::vector2_aligned<Tp, Alloc>& y)
  {
    x.swap(x);
  }
};
#endif
