// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__BICHART__HPP__
#define __UTILS__BICHART__HPP__ 1

//
// a compact implementation of bi-chart structure...
//

#include <vector>

namespace utils
{
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class bichart
  {
  public:
    typedef Tp value_type;
    typedef std::vector<Tp, Alloc> chart_type;
    
    typedef typename chart_type::size_type        size_type;
    typedef typename chart_type::difference_type  difference_type;
    
    typedef typename chart_type::const_reference         const_reference;
    typedef typename chart_type::reference               reference;
    
    typedef typename chart_type::pointer                 pointer;
    
  public:
    bichart() : __chart(), __size1(0), __size2(0) {}
    bichart(size_type s1, size_type s2) : __chart(__chart_size(s1, s2)), __size1(s1), __size2(s2) {}
    bichart(size_type s1, size_type s2, const value_type& v) : __chart(__chart_size(s1, s2), v), __size1(s1), __size2(s2) {}
    
  public:
    template <typename _Tp, typename _Alloc>
    void swap(bichart<_Tp, _Alloc>& x)
    {
      __chart.swap(x.__chart);
      std::swap(__size1, x.__size1);
      std::swap(__size2, x.__size2);
    }
    
    bool empty() const { return __chart.empty(); }
    size_type size1() const { return __size1; }
    size_type size2() const { return __size2; }
    
    void clear() { __chart.clear(); __size1 = 0; __size2 = 0; }
    void resize(size_type s1, size_type s2) { __chart.resize(__chart_size(s1, s2)); __size1 = s1; __size2 = s2; }
    void resize(size_type s1, size_type s2, const value_type& v) { __chart.resize(__chart_size(s1, s2), v); __size1 = s1; __size2 = s2; }
    
    void reserve(size_type s1, size_type s2) { __chart.reserve(__chart_size(s1, s2)); }
    
    inline const_reference operator()(size_type x1, size_type x2, size_type x3, size_type x4) const
    {
      return __chart[__position(x1, x2, x3, x4, __size1, __size2)];
    }
    
    inline reference       operator()(size_type x1, size_type x2, size_type x3, size_type x4)
    {
      return __chart[__position(x1, x2, x3, x4, __size1, __size2)];
    }
    
  private:
    static inline
    size_type __chart_size(size_type s)
    {
      return s * ((s >> 1) + 1);
    }
    
    static inline
    size_type __chart_size(size_type s1, size_type s2)
    {
      return __chart_size(s1) * __chart_size(s2);
    }

    static inline
    size_type __position(size_type x1, size_type x2, size_type x3, size_type x4, size_type s1, size_type s2)
    {
      return __position(x2 - x1, x2, s1) * __chart_size(s2) + __position(x4 - x3, x4, s2);
    }
    
    static inline
    size_type __position(size_type range, size_type pos, size_type size)
    {
      const size_type offset = (size >> 1) + 1;
      const size_type offset_mask = (range >= offset) - 1;
      const size_type offset_inv = ~offset_mask;
      
      const size_type x1_pos = (offset_inv & (size - range)) | (offset_mask & range);
      const size_type x2_pos = (offset_inv & (pos - range))  | (offset_mask & pos);
      
      return x1_pos * size + x2_pos;
    }
    
  private:
    chart_type __chart;
    size_type  __size1;
    size_type  __size2;
  };
};

namespace std
{
  template <typename _Tp, typename _Alloc>
  inline
  void swap(utils::bichart<_Tp, _Alloc>& x, utils::bichart<_Tp, _Alloc>& y)
  {
    x.swap(y);
  }
  
};


#endif
