// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__ALLOCINFO_ALLOCATOR__HPP__
#define __UTILS__ALLOCINFO_ALLOCATOR__HPP__ 1

#include <stdint.h>

#include <boost/thread.hpp>

#include <utils/config.hpp>

namespace utils
{
  struct allocinfo
  {
    typedef int64_t allocated_type;
    
    allocinfo() {}
    
    allocated_type& allocated()
    {
#ifdef HAVE_TLS
      static __thread allocated_type __allocated = 0;
      return __allocated;
#else
      static boost::thread_specific_ptr<allocated_type> __allocated;
      
      if (! __allocated.get())
	__allocated.reset(new allocated_type(0));
      return *__allocated;
#endif
    }

  };
  
  template <typename _Tp, typename _Alloc=std::allocator<_Tp> >
  class allocinfo_allocator
  {
  public:
    typedef size_t     size_type;
    typedef ptrdiff_t  difference_type;
    typedef _Tp*       pointer;
    typedef const _Tp* const_pointer;
    typedef _Tp&       reference;
    typedef const _Tp& const_reference;
    typedef _Tp        value_type;
    typedef allocinfo::allocated_type   allocated_type;
    
    template<typename _Tp1>
    struct rebind
    { typedef allocinfo_allocator<_Tp1, _Alloc> other; };
    
  private:
    typedef typename _Alloc::template rebind<_Tp>::other allocator_type;
    typedef allocinfo allocinfo_type;
    
  public:
    
    allocinfo_allocator() throw() {}
    allocinfo_allocator(const allocinfo_allocator& ) throw() {}
    template <typename _Tp1>
    allocinfo_allocator(const allocinfo_allocator<_Tp1, _Alloc>&) throw() {}
    template <typename _Tp1, typename _Alloc1>
    allocinfo_allocator(const allocinfo_allocator<_Tp1, _Alloc1>&) throw() {}
    ~allocinfo_allocator() throw() {}
    
    pointer address(reference __x) const { return &__x; }
    const_pointer address(const_reference __x) const { return &__x; }
    
    
    pointer allocate(size_type __n, const void* __p = 0)
    {
      allocated() += sizeof(_Tp) * __n;
      return allocator().allocate(__n, __p);
    }
    
    void deallocate(pointer __p, size_type __n)
    { 
      allocated() -= sizeof(_Tp) * __n;
      allocator().deallocate(__p, __n);
    }
    
    size_type max_size() const throw()
    { return size_t(-1) / sizeof(_Tp); }
    

    void construct(pointer __p, const _Tp& __val)
    { ::new(__p) _Tp(__val); }
    
    
    void destroy(pointer __p) { __p->~_Tp(); }

  private:
    allocinfo_type& _allocinfo()
    {
      static allocinfo_type __allocinfo;
      return __allocinfo;
    }
    
    allocated_type& allocated()
    {
      return _allocinfo().allocated();
    }
    
    allocator_type& allocator()
    {
      static allocator_type __allocator;
      return __allocator;
    }
  };
  
  template<typename _Tp, typename _Alloc>
  inline bool
  operator==(const allocinfo_allocator<_Tp, _Alloc>&, const allocinfo_allocator<_Tp, _Alloc>&)
  { return true; }
  
  template<typename _Tp, typename _Alloc>
  inline bool
  operator!=(const allocinfo_allocator<_Tp, _Alloc>&, const allocinfo_allocator<_Tp, _Alloc>&)
  { return false; }
  
};

#endif
