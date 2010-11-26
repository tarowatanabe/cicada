// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__VECTOR3__HPP__
#define __UTILS__VECTOR3__HPP__ 1

#include <vector>

namespace utils
{
  
  template <typename _Tp, typename _Alloc=std::allocator<_Tp> >
  class vector3
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
    
  public:
    vector3()
      : __base(), __size1(), __size2(), __size3() {}
    vector3(size_type __s1, size_type __s2, size_type __s3)
      : __base(__s1 * __s2 * __s3), __size1(__s1), __size2(__s2), __size3(__s3) {}
    vector3(size_type __s1, size_type __s2, size_type __s3, const value_type& value)
      : __base(__s1 * __s2 * __s3, value), __size1(__s1), __size2(__s2), __size3(__s3) {}
    
    bool empty() const { return __base.empty(); }
    size_type size1() const { return __size1; }
    size_type size2() const { return __size2; }
    size_type size3() const { return __size3; }
    
    void reserve(size_type __s1, size_type __s2, size_type __s3)
    {
      __base.reserve(__s1 * __s2 * __s3);
    }
    
    void resize(size_type __s1, size_type __s2, size_type __s3)
    {
      __size1 = __s1;
      __size2 = __s2;
      __size3 = __s3;
      __base.resize(__s1 * __s2 * __s3);
    }
    void resize(size_type __s1, size_type __s2, size_type __s3, const value_type& value)
    {
      __size1 = __s1;
      __size2 = __s2;
      __size3 = __s3;
      __base.resize(__s1 * __s2 * __s3, value);
    }
    void clear()
    {
      __base.clear();
      __size1 = 0;
      __size2 = 0;
      __size3 = 0;
    }
    void swap(vector3& x)
    {
      using namespace std;
      __base.swap(x.__base);
      swap(__size1, x.__size1);
      swap(__size2, x.__size2);
      swap(__size3, x.__size3);
    }
    
    inline const_reference operator()(size_type pos1, size_type pos2, size_type pos3) const
    { 
      return __base[pos1 * __size2 * __size3 + pos2 * __size3 + pos3];
    }
    inline       reference operator()(size_type pos1, size_type pos2, size_type pos3)
    {
      return __base[pos1 * __size2 * __size3 + pos2 * __size3 + pos3];
    }
    
    inline const_iterator begin() const { return __base.begin(); }
    inline       iterator begin()       { return __base.begin(); }
    
    inline const_iterator end() const { return __base.end(); }
    inline       iterator end()       { return __base.end(); }
    
    inline const_iterator begin(size_type pos) const { return __base.begin() + pos * __size2 * __size3; }
    inline       iterator begin(size_type pos)       { return __base.begin() + pos * __size2 * __size3; }
    
    inline const_iterator end(size_type pos) const { return __base.begin() + (pos + 1) * __size2 * __size3; }
    inline       iterator end(size_type pos)       { return __base.begin() + (pos + 1) * __size2 * __size3; }
    
    inline const_iterator begin(size_type pos1, size_type pos2) const { return __base.begin() + pos1 * __size2 * __size3 + pos2 * __size3; }
    inline       iterator begin(size_type pos1, size_type pos2)       { return __base.begin() + pos1 * __size2 * __size3 + pos2 * __size3; }
    
    inline const_iterator end(size_type pos1, size_type pos2) const { return __base.begin() + pos1 * __size2 * __size3 + (pos2 + 1) * __size3; }
    inline       iterator end(size_type pos1, size_type pos2)       { return __base.begin() + pos1 * __size2 * __size3 + (pos2 + 1) * __size3; }

  private:
    base_type __base;
    size_type __size1;
    size_type __size2;
    size_type __size3;
  };
  
};

namespace std
{
  template <typename Tp, typename Alloc>
  inline
  void swap(utils::vector3<Tp, Alloc>& x, utils::vector3<Tp, Alloc>& y)
  {
    x.swap(x);
  }
};
#endif
