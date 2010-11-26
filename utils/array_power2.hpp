// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__ARRAY_POWER2__HPP__
#define __UTILS__ARRAY_POWER2__HPP__ 1

#include <algorithm>
#include <memory>
#include <utility>
#include <stdexcept>
#include <iterator>

#include <utils/bithack.hpp>
#include <utils/memory.hpp>

//
// it is an implementation of array, but performs no initialization..
//

namespace utils
{
  
  template <typename _Tp, size_t Size, typename _Alloc>
  struct __array_power2_base : public _Alloc
  {
    typedef _Tp        value_type;
    typedef size_t     size_type;
    typedef _Tp*       pointer;
    
    typedef _Tp*       iterator;
    typedef const _Tp* const_iterator;
    
    typedef _Tp&       reference;
    typedef const _Tp& const_reference;
    
    typedef std::reverse_iterator<iterator>       reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    
    pointer m_base;
    
    __array_power2_base() : m_base() { allocate(); }
    ~__array_power2_base() { deallocate(); }
    
    void allocate()
    {
      m_base = static_cast<_Alloc&>(*this).allocate(Size);
    }
    
    void deallocate()
    {
      if (m_base)
	static_cast<_Alloc&>(*this).deallocate(m_base, Size);
    }
  };
  
  template <typename _Tp, size_t _Size, typename _Alloc=std::allocator<_Tp> >
  class array_power2
  {
  private:
    static const bool   __is_power2 = utils::bithack::static_is_power2<_Size>::result;
    static const size_t __next_power2 = size_t(utils::bithack::static_next_largest_power2<_Size>::result);
    static const size_t static_size = (_Size == 0 ? size_t(1) : (__is_power2 ? _Size : __next_power2));
    
  private:
    typedef __array_power2_base<_Tp, static_size, _Alloc> base_type;
    
  public:
    typedef typename base_type::value_type    value_type;
    typedef typename base_type::size_type     size_type;

    typedef typename base_type::iterator               iterator;
    typedef typename base_type::const_iterator         const_iterator;
    typedef typename base_type::reverse_iterator       reverse_iterator;
    typedef typename base_type::const_reverse_iterator const_reverse_iterator;    
    typedef typename base_type::reference              reference;
    typedef typename base_type::const_reference        const_reference;
    
    typedef array_power2<_Tp,_Size,_Alloc> self_type;
    
  public:
    array_power2() { std::uninitialized_fill(begin(), end(), value_type()); }
    array_power2(const array_power2& x) { std::uninitialized_copy(x.begin(), x.end(), begin()); }
    ~array_power2() { utils::destroy_range(begin(), end()); }
    
    array_power2& operator=(const array_power2& x)
    {
      std::copy(x.begin(), x.end(), begin());
      return *this;
    }
    
    void assign(const self_type& x)
    {
      std::copy(x.begin(), x.end(), begin());
    }
    
  public:
    bool empty() const { return false; }
    const size_type size() const { return static_size; }
    void clear() { std::fill(begin(), end(), value_type()); }
    
    iterator begin() { return base.m_base; }
    const_iterator begin() const { return base.m_base; }
    iterator end() { return base.m_base + static_size; }
    const_iterator end() const { return base.m_base + static_size; }
    
    reverse_iterator rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const { return reverse_iterator(end()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const { return reverse_iterator(begin()); }

    reference front() { return *(base.m_base); }
    const_reference front() const { return *(base.m_base); }
    
    reference back() { return *(base.m_base + static_size - 1); }
    const_reference back() const { return *(base.m_base + static_size - 1); }
    
    reference operator[](size_type x) { return base.m_base[x]; }
    const_reference operator[](size_type x) const { return base.m_base[x]; }
    
    void swap(self_type& x)
    {
      std::swap_ranges(begin(), end(), x.begin());
    }

  private:
    
    
  private:    
    base_type base;
  };

  template <typename Tp, size_t Size, typename Alloc>
  inline
  bool operator<(const array_power2<Tp,Size,Alloc>& x,
		 const array_power2<Tp,Size,Alloc>& y)
  {
    return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
  }
  
  template <typename Tp, size_t Size, typename Alloc>
  inline
  bool operator>(const array_power2<Tp,Size,Alloc>& x,
		 const array_power2<Tp,Size,Alloc>& y)
  {
    return y < x;
  }

  template <typename Tp, size_t Size, typename Alloc>
  inline
  bool operator<=(const array_power2<Tp,Size,Alloc>& x,
		  const array_power2<Tp,Size,Alloc>& y)
  {
    return ! (y < x);
  }
  
  template <typename Tp, size_t Size, typename Alloc>
  inline
  bool operator>=(const array_power2<Tp,Size,Alloc>& x,
		  const array_power2<Tp,Size,Alloc>& y)
  {
    return ! (x < y);
  }
  
  template <typename Tp, size_t Size, typename Alloc>
  inline
  bool operator==(const array_power2<Tp,Size,Alloc>& x,
		  const array_power2<Tp,Size,Alloc>& y)
  {
    return std::equal(x.begin(), x.end(), y.begin());
  }
  
  template <typename Tp, size_t Size, typename Alloc>
  inline
  bool operator!=(const array_power2<Tp,Size,Alloc>& x,
		  const array_power2<Tp,Size,Alloc>& y)
  {
    return ! (x == y);
  }
  
};

namespace std
{
  template <typename Tp, size_t Size, typename Alloc>
  inline
  void swap(const utils::array_power2<Tp,Size,Alloc>& x,
	    const utils::array_power2<Tp,Size,Alloc>& y)
  {
    x.swap(y);
  }
};

#endif
