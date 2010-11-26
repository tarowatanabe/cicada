// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__ALIGNED_ALLOCATOR__HPP__
#define __UTILS__ALIGNED_ALLOCATOR__HPP__ 1

// return 16-byte aligned memory region...
// by default, malloc will be suffice in most machines.

#include <stdexcept>

#include <string>
#include <algorithm>

#include <cstdlib>
#include <cstring>

#include <utils/config.hpp>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#ifdef HAVE_MM_MALLOC_H
#include <mm_malloc.h>
#endif


namespace utils
{
  
  // alignment must be power of two
  
  template <typename _Tp>
  struct aligned_allocator
  {
    
  public:
    typedef size_t     size_type;
    typedef ptrdiff_t  difference_type;
    typedef _Tp*       pointer;
    typedef const _Tp* const_pointer;
    typedef _Tp&       reference;
    typedef const _Tp& const_reference;
    typedef _Tp        value_type;
    
    template<typename _Tp1>
    struct rebind
    { typedef aligned_allocator<_Tp1> other; };
    
  public:
    aligned_allocator() throw() { }
    aligned_allocator(const aligned_allocator&) throw() { }
    template<typename _Tp1>
    aligned_allocator(const aligned_allocator<_Tp1>&) throw() { }
    ~aligned_allocator() throw() { }
    
  public:
    pointer
    address(reference __x) const { return &__x; }
    
    const_pointer
    address(const_reference __x) const { return &__x; }
    
    size_type
    max_size() const throw() 
    { return size_type(-1) / sizeof(_Tp); }
    
    void 
    construct(pointer __p, const _Tp& __val) 
    { ::new(__p) _Tp(__val); }
    
    void 
    destroy(pointer __p) { __p->~_Tp(); }
    
    pointer
    allocate(size_type __n, const void* tmptmp= 0)
    {
#if defined(HAVE_POSIX_MEMALIGN)
      void* p = 0;
      const int err = ::posix_memalign(&p, 16, sizeof(_Tp) * __n);
      if (err != 0 || p == 0)
	throw std::bad_alloc();
      if (size_type(p) & size_type(16 - 1) != 0)
	throw std::runtime_error("not aligned");
      return static_cast<pointer>(p);
#elif defined(HAVE_MEMALIGN)
      // I'm not sure whether the memaligned block can be freeed... in glibc, it can
      // but not specified in Sun's man page
      void* p = ::memalign(16, sizeof(_Tp) * __n);
      if (p == 0)
	throw std::bad_alloc();
      if (size_type(p) & size_type(16 - 1) != 0)
	throw std::runtime_error("not aligned");
      return static_cast<pointer>(p);
#elif defined(HAVE_MM_MALLOC_H)
#define __UTILS__ALIGNED_ALLOCATOR__MALLOC_BY_MM_ALLOC__ 1
      void* p = _mm_malloc(sizeof(_Tp) * __n, 16);
      if (p == 0)
	throw std::bad_alloc();
      return static_cast<pointer>(p);
#else
      void* p = ::malloc(sizeof(_Tp) * __n);
      if (p == 0)
	throw std::bad_alloc();
      if (size_type(p) & size_type(16 - 1) != 0)
	throw std::runtime_error("not aligned");
      return static_cast<pointer>(p);
#endif
    }
    
    void
    deallocate(pointer __p, size_type __n)
    {
#ifdef __UTILS__ALIGNED_ALLOCATOR__MALLOC_BY_MM_ALLOC__
      if (__p)
	_mm_free(__p);
#else
      if (__p)
	::free(__p);
#endif
    }
  };

  template <typename _Tp>
  inline bool operator==(const aligned_allocator<_Tp>& x,
                         const aligned_allocator<_Tp>& y)
  { return true; }
  
  template <typename _Tp>
  inline bool operator!=(const aligned_allocator<_Tp>& x,
                         const aligned_allocator<_Tp>& y)
  { return false; }

};

#endif
