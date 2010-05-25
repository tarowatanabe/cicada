// -*- mode: c++ -*-

#ifndef __UTILS__MPI_ALLOCATOR__HPP__
#define __UTILS__MPI_ALLOCATOR__HPP__ 1

#include <new>
#include <stdexcept>
#include <string>
#include <algorithm>

#include <mpi.h>

namespace utils
{
  template <typename _Tp>
  struct mpi_allocator
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
    { typedef mpi_allocator<_Tp1> other; };
    
  public:
    mpi_allocator() throw() { }
    mpi_allocator(const mpi_allocator&) throw() { }
    
    template<typename _Tp1>
    mpi_allocator(const mpi_allocator<_Tp1>&) throw() { }
    
    ~mpi_allocator() throw() { }
  public:
    pointer
    address(reference __x) const { return &__x; }
    
    const_pointer
    address(const_reference __x) const { return &__x; }
    
    size_type
    max_size() const throw() 
    { return size_t(-1) / sizeof(_Tp); }
    
    void 
    construct(pointer __p, const _Tp& __val) 
    { ::new(__p) _Tp(__val); }
    
    void 
    destroy(pointer __p) { __p->~_Tp(); }
    
    pointer
    allocate(size_type __n, const void* tmptmp= 0)
    {
      if (__n == 0) return 0;
      pointer p = static_cast<pointer>(MPI::Alloc_mem(sizeof(_Tp) * __n, MPI_INFO_NULL));
      if (p == 0)
	throw std::bad_alloc();
      return p;
    }
    
    void
    deallocate(pointer __p, size_type __n)
    {
      if (__p)
	MPI::Free_mem(__p);
    }
  };
};

#endif
