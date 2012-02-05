// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MULVECTOR2__HPP__
#define __UTILS__MULVECTOR2__HPP__ 1

#include <vector>
#include <boost/swap.hpp>

namespace utils
{
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class mulvector2
  {
  public:
    typedef size_t      size_type;
    typedef ptrdiff_t   difference_type;
    
  private:
    typedef typename Alloc::template rebind<size_type>::other offset_alloc_type;
    
    typedef std::vector<Tp, Alloc>                    impl_type;
    typedef std::vector<size_type, offset_alloc_type> offset_set_type;
    
  private:
    struct __reference
    {
    public:
      typedef typename impl_type::size_type       size_type;
      typedef typename impl_type::difference_type difference_type;
      typedef typename impl_type::value_type      value_type;
      
      typedef typename impl_type::const_iterator          const_iterator;
      typedef typename impl_type::const_reverse_iterator  const_reverse_iterator;
      typedef typename impl_type::const_reference         const_reference;
      
    public:
      __reference(const_iterator __first, const_iterator __last) : first(__first), last(__last) {}
      
    public:
      const_reference operator[](size_type x) const{ return *(first + x); }
      const_reference front() const { return *first; }
      const_reference back() const { return *(last - 1); }
      
      const_iterator begin() const{ return first; }
      const_iterator end() const { return last; }
      
      size_type size() const { return std::distance(first, last); }
      bool empty() const { return first == last; }
      
    public:
      friend bool operator==(const __reference& x, const __reference& y) { return x.first == y.first && x.last == y.last; }
      friend bool operator!=(const __reference& x, const __reference& y) { return x.first != y.first || x.last != y.last; }
      friend bool operator<(const __reference& x, const __reference& y) { return x.first < y.first || (!(y.first < x.first) && x.last < y.last); } 
      friend bool operator>(const __reference& x, const __reference& y) { return y < x; } 
      friend bool operator<=(const __reference& x, const __reference& y) { return ! (y < x); }
      friend bool operator>=(const __reference& x, const __reference& y) { return ! (x < y); }
      
    private:
      const_iterator first;
      const_iterator last;
    };

    struct __iterator
    {
    public:
      typedef typename offset_set_type::const_iterator iter_type;
      typedef __reference const_reference;

    public:
      typedef typename offset_set_type::size_type       size_type;
      typedef typename offset_set_type::difference_type difference_type;
      
    public:
      __iterator(iter_type __iter, const impl_type& __impl) : iter(__iter), pimpl(&__impl) {}
      
      const_reference operator*() const
      {
	return const_reference(pimpl->begin() + *(iter - 1), pimpl->begin() + *iter);
      }
      
      const_reference operator[](difference_type __n) const
      {
	return *(*this + __n);
      }
      
      
      __iterator& operator++() { ++ iter; return *this; }
      __iterator& operator--() { -- iter; return *this; }
      
      __iterator& operator+=(difference_type __n) { iter += __n; return *this; }
      __iterator& operator-=(difference_type __n) { iter -= __n; return *this; }

      __iterator operator++(int)
      {
	__iterator __tmp = *this;
	++ *this;
	return __tmp;
      }
      
      __iterator operator--(int)
      {
	__iterator __tmp = *this;
	-- *this;
	return __tmp;
      }

      __iterator operator+(difference_type __n) const
      {
	__iterator __tmp = *this;
	return __tmp += __n;
      }

      __iterator operator-(difference_type __n) const
      {
	__iterator __tmp = *this;
	return __tmp -= __n;
      }
      
      friend bool operator==(const __iterator& x, const __iterator& y) { return x.iter == y.iter && x.pimpl == y.pimpl; }
      friend bool operator!=(const __iterator& x, const __iterator& y) { return x.iter != y.iter || x.pimpl != y.pimpl; }
      friend bool operator<(const __iterator& x, const __iterator& y) { return x.pimpl < y.pimpl || (x.pimpl == y.pimpl && x.iter < y.iter); } 
      friend bool operator>(const __iterator& x, const __iterator& y) { return y < x; } 
      friend bool operator<=(const __iterator& x, const __iterator& y) { return ! (y < x); }
      friend bool operator>=(const __iterator& x, const __iterator& y) { return ! (x < y); }
      
    private:
      iter_type         iter;
      const impl_type*  pimpl;
    };

  public:
    typedef __reference const_reference;
    typedef __iterator  const_iterator;
    
  public:
    mulvector2() : impl(), offsets(1, 0) {}

  public:
    const_reference operator[](size_type x) const { return const_reference(impl.begin() + offsets[x], impl.begin() + offsets[x + 1]); }
    const_reference front() const { return const_reference(impl.begin() + offsets[0], impl.begin() + offsets[1]); }
    const_reference back() const { return const_reference(impl.begin() + offsets[offsets.size() - 2], impl.begin() + offsets[offsets.size() - 1]); }
    
    const_iterator begin() const { return const_iterator(offsets.begin() + 1, impl); }
    const_iterator end() const { return const_iterator(offsets.end(), impl); }
    
    template <typename Iterator>
    size_type push_back(Iterator first, Iterator last)
    {
      impl.insert(impl.end(), first, last);
      offsets.push_back(impl.size());
      
      return size() - 1;
    }
    
    size_type push_back()
    {
      offsets.push_back(impl.size());
      return size() - 1;
    }
    
    void reserve(size_type size1, size_type size2)
    {
      offsets.reserve(size1 + 1);
      impl.reserve(size1 * size2);
    }
    
    void clear()
    {
      impl.clear();
      offsets.resize(1);
    }
    
    size_type size() const { return offsets.size() - 1; }
    bool empty() const { return offsets.size() == 1; }

    void swap(mulvector2& x)
    {
      impl.swap(x.impl);
      offsets.swap(x.offsets);
    }
    
  private:
    impl_type       impl;
    offset_set_type offsets;
  };
  
};

namespace std
{
  template <typename Tp, typename Alloc>
  inline
  void swap(utils::mulvector2<Tp,Alloc>& x, utils::mulvector2<Tp,Alloc>& y)
  {
    x.swap(y);
  }
};

#endif
