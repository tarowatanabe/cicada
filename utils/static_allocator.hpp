// -*- mode: c++ -*-

#ifndef __UTILS__STATIC_ALLOCATOR__HPP__
#define __UTILS__STATIC_ALLOCATOR__HPP__ 1

#include <stdexcept>
#include <climits>


#include <utils/config.hpp>
#include <boost/thread.hpp>

namespace utils
{
  template <size_t _Size, typename _Alloc, bool Bool >
  struct __static_allocator_impl
  {
    
  };
  
  template <size_t _Size, typename _Alloc >
  struct __static_allocator_impl<_Size, _Alloc, true> 
  {
  public:
    typedef size_t      size_type;
    typedef ptrdiff_t   difference_type;
    typedef char*       pointer;
    typedef const char* const_pointer;
    typedef char&       reference;
    typedef const char& const_reference;
    typedef char        value_type;
    
    typedef __static_allocator_impl<_Size, _Alloc, false> self_type;

    typedef typename _Alloc::template rebind<char>::other allocator_type;
    
    struct Pool
    {
      pointer pool;
      short size;
      Pool() : pool(), size(0) {}
      ~Pool() {
	while (pool) {
	  pointer p = pool;
	  pool = *((pointer*) p);
	  alloc().deallocate(p, _Size);
	}
      }
    };
    typedef Pool pool_type;
    
  public:
    __static_allocator_impl() {}
    __static_allocator_impl(const self_type& x) {}
    self_type& operator=(const self_type& x)
    {
      return *this;
    }
    ~__static_allocator_impl() throw() { }
    
    pointer allocate()
    {
      // access to local pool
      pool_type& pool = local_pool();
      if (pool.pool) {
	pointer p = pool.pool;
	pool.pool = *((pointer*) p);
	-- pool.size;
	return p;
      } else
	return alloc().allocate(_Size);
    }
    
    void deallocate(pointer p)
    {
      if (! p) return;
      
      pool_type& pool = local_pool();
      if (pool.size < 256) {
	*((pointer*) p) = pool.pool;
	pool.pool = p;
	++ pool.size;
      } else
	alloc().deallocate(p, _Size);
    }
    
  private:
    
    pool_type& local_pool()
    {
    
      // differentiate here for __thread keyword... is this correct?
#ifdef HAVE_TLS
      static __thread pool_type* local_pool = 0;
      static boost::thread_specific_ptr<pool_type> local_pool_tss;
      
      if (local_pool == 0) {
	local_pool_tss.reset(new pool_type());
	local_pool = local_pool_tss.get();
      }
      return *local_pool;
#else
      static boost::thread_specific_ptr<pool_type> local_pool;
      if (! local_pool.get())
	local_pool.reset(new pool_type());
      return *local_pool;
#endif
    }
    
    static allocator_type& alloc()
    {
      static allocator_type __alloc;
      return __alloc;
    }    
  };
  

  // we do not do any recycling...
  template <size_t _Size, typename _Alloc >
  struct __static_allocator_impl<_Size, _Alloc, false>
  {
  public:
    typedef size_t      size_type;
    typedef ptrdiff_t   difference_type;
    typedef char*       pointer;
    typedef const char* const_pointer;
    typedef char&       reference;
    typedef const char& const_reference;
    typedef char        value_type;
    
    typedef typename _Alloc::template rebind<char>::other allocator_type;
    typedef __static_allocator_impl<_Size, _Alloc, false> self_type;
    
  public:
    __static_allocator_impl() {}
    __static_allocator_impl(const self_type& x) {}
    self_type& operator=(const self_type& x)
    {
      return *this;
    }
    ~__static_allocator_impl() throw() { }
    
    pointer allocate()
    {
      return alloc().allocate(_Size);
    }
    
    void deallocate(pointer p)
    {
      if (! p) return;
      
      alloc().deallocate(p, _Size);
    }
    
  private:
    
    static allocator_type& alloc()
    {
      static allocator_type __alloc;
      return __alloc;
    }
  };


  template <size_t _Size, typename _Alloc >
  struct __static_allocator_base
  {
    // use of static_allocator if the size is larger than 128
    // we set 256, that is the chunk-size for mt-allocator or tcmalloc
    typedef __static_allocator_impl<_Size, _Alloc, _Size >= sizeof(void*) > base_type;
    
    base_type& base()
    {
      static base_type __base;
      return __base;
    }
  };

  template <typename _Tp, size_t _Size, typename _Alloc=std::allocator<_Tp> >
  struct static_allocator : public __static_allocator_base<sizeof(_Tp) * _Size, _Alloc>
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
    { typedef static_allocator<_Tp1, _Size, _Alloc> other; };
    
    using __static_allocator_base<sizeof(_Tp) * _Size, _Alloc>::base;

  private:
    typedef typename _Alloc::template rebind<_Tp >::other base_allocator_type;
    
  public:
    static_allocator() throw() { }
    
    static_allocator(const static_allocator&) throw() { }
    
    template<typename _Tp1>
    static_allocator(const static_allocator<_Tp1,_Size,_Alloc>&) throw() { }
    
    ~static_allocator() throw() { }
    
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
    allocate(size_type __n, const void* = 0)
    {
      if (__n != _Size)
	return base_allocator().allocate(__n);
      else
	return (pointer) base().allocate();
    }
    
    void
    deallocate(pointer __p, size_type __n)
    {
      if (__n != _Size)
	base_allocator().deallocate(__p, __n);
      else
	base().deallocate((char*) __p);
    }
    
  private:
    base_allocator_type& base_allocator()
    {
      static base_allocator_type _alloc;
      return _alloc;
    }

  };

  template <typename _Tp, size_t _Size, typename _Alloc>
  inline bool operator==(const static_allocator<_Tp,_Size,_Alloc>& x,
                         const static_allocator<_Tp,_Size,_Alloc>& y)
  { return true; }
  
  template <typename _Tp, size_t _Size, typename _Alloc>
  inline bool operator!=(const static_allocator<_Tp,_Size,_Alloc>& x,
                         const static_allocator<_Tp,_Size,_Alloc>& y)
  { return false; }

};

#endif
