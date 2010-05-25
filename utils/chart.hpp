// -*- mode: c++ -*-

#ifndef __UTILS__CHART__HPP__
#define __UTILS__CHART__HPP__ 1

//
// a compact implementation of chart structure...
//

#include <vector>

namespace utils
{
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class chart
  {
  public:
    typedef Tp value_type;
    typedef std::vector<Tp, Alloc> chart_type;
    
    typedef typename chart_type::size_type        size_type;
    typedef typename chart_type::difference_type  difference_type;

    typedef typename chart_type::const_reference         const_reference;
    typedef typename chart_type::reference               reference;
    
    typedef typename chart_type::const_iterator          const_iterator;
    typedef typename chart_type::iterator                iterator;
    
    typedef typename chart_type::const_reverse_iterator  const_reverse_iterator;
    typedef typename chart_type::reverse_iterator        reverse_iterator;
    
    typedef typename chart_type::pointer                 pointer;
  public:
    chart() : __chart(), __size(0) {}
    chart(size_type __s) : __chart(__chart_size(__s)), __size(__s) {}
    chart(size_type __s, const value_type& __v) : __chart(__chart_size(__s), __v), __size(__s) {}
    
  public:
    template <typename _Tp, typename _Alloc>
    void swap(chart<_Tp, _Alloc>& x) {
      __chart.swap(x.__chart);
      std::swap(__size, x.__size);
    }
    
    bool empty() const { return __chart.empty(); }
    size_type size() const { return __size; }
    
    void clear() { __chart.clear(); __size = 0; }
    void resize(size_type __s) { __chart.resize(__chart_size(__s)); __size = __s; }
    void resize(size_type __s, const value_type& __v) { __chart.resize(__chart_size(__s), __v); __size = __s; }
    void reserve(size_type __s) { __chart.reserve(__chart_size(__s)); }
    
    inline const_reference operator()(size_type x1, size_type x2) const { return __chart[__position(x2 - x1, x2)]; }
    inline reference       operator()(size_type x1, size_type x2)       { return __chart[__position(x2 - x1, x2)]; }
    
    inline const_iterator begin(size_type length) const { return __chart.begin() + __position(length, length); }
    inline iterator       begin(size_type length)       { return __chart.begin() + __position(length, length); }
    
    inline const_iterator end(size_type length) const { return __chart.begin() + __position(length, __size); }
    inline iterator       end(size_type length)       { return __chart.begin() + __position(length, __size); }
    
    
    inline const_iterator begin() const { return __chart.begin(); }
    inline iterator       begin()       { return __chart.begin(); }
    
    inline const_iterator end() const { return __chart.end(); }
    inline iterator       end()       { return __chart.end(); }
    
    inline const_reverse_iterator rbegin() const { return __chart.rbegin(); }
    inline reverse_iterator       rbegin()       { return __chart.rbegin(); }
    
    inline const_reverse_iterator rend() const { return __chart.rend(); }
    inline reverse_iterator       rend()       { return __chart.rend(); }
    
  private:
    static inline
    size_type __chart_size(size_type __s)
    {
      return __s * ((__s >> 1) + 1);
    }
    
    inline
    size_type __position(size_type range, size_type pos) const
    {
      const size_type offset = (__size >> 1) + 1;
      const size_type offset_mask = (range >= offset) - 1;
      const size_type offset_inv = ~offset_mask;
      
      const size_type x1_pos = (offset_inv & (__size - range)) | (offset_mask & range);
      const size_type x2_pos = (offset_inv & (pos - range))    | (offset_mask & pos);
      
      return x1_pos * __size + x2_pos;
    }
    
  private:
    chart_type __chart;
    size_type  __size;
  };
};

namespace std
{
  template <typename _Tp, typename _Alloc>
  inline
  void swap(utils::chart<_Tp, _Alloc>& x, utils::chart<_Tp, _Alloc>& y)
  {
    x.swap(y);
  }
  
};


#endif
