// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SNAPPY_VECTOR__HPP__
#define __UTILS__SNAPPY_VECTOR__HPP__ 1

// snappy vector with random access support

#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>

#include <boost/filesystem.hpp>
#include <boost/array.hpp>

#include <utils/repository.hpp>
#include <utils/filesystem.hpp>
#include <utils/bithack.hpp>
#include <utils/spinlock.hpp>

#include <utils/snappy_file.hpp>

namespace utils
{
  template <typename Tp>
  struct __snappy_vector_base
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef char      char_type;
    typedef char      byte_type;
    typedef uint64_t  off_type;
    
    typedef Tp value_type;
    
    typedef boost::filesystem::path path_type;
    
    static const size_type __buffer_non_aligned_size = (128 / sizeof(value_type) <= 1 ? 1 : 128 / sizeof(value_type));
    static const size_type __buffer_non_aligned_is_power2 = utils::bithack::static_is_power2<__buffer_non_aligned_size>::value;
    static const size_type __buffer_non_aligned_next_power2 = utils::bithack::static_next_largest_power2<__buffer_non_aligned_size>::value;
    
    static const size_type __buffer_size = (__buffer_non_aligned_is_power2 ? __buffer_non_aligned_size : __buffer_non_aligned_next_power2);
    static const size_type __buffer_mask = __buffer_size - 1;
    static const size_type __buffer_shift = utils::bithack::static_bit_count<__buffer_mask>::value;
  };
  

  template <typename Tp, typename Impl>
  struct __snappy_vector_iterator
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Tp        value_type;
    typedef Impl      impl_type;

    typedef __snappy_vector_iterator<Tp,Impl> self_type;
    
    typedef const Tp* pointer;
    typedef const Tp& reference;
    typedef std::random_access_iterator_tag   iterator_category;
    
    __snappy_vector_iterator(size_type pos, const impl_type& impl)
      : __pos(pos), __impl(&impl) {}
    __snappy_vector_iterator() : __pos(0), __impl(0) {}
    
    const value_type& operator*() const { return __impl->operator[](__pos); }
    const value_type* operator->() const { return &__impl->operator[](__pos); }
    
    self_type& operator++() { ++ __pos; return *this; }
    self_type& operator--() { -- __pos; return *this; }
    
    self_type& operator+=(difference_type __n) { __pos += __n; return *this; }
    self_type& operator-=(difference_type __n) { __pos -= __n; return *this; }
    
    self_type operator++(int) { self_type __tmp = *this; ++ *this; return __tmp; }
    self_type operator--(int) { self_type __tmp = *this; -- *this; return __tmp; }
    self_type operator+(difference_type __n) const { self_type __tmp = *this; return __tmp += __n; }
    self_type operator-(difference_type __n) const { self_type __tmp = *this; return __tmp -= __n; }
    
    size_type        __pos;
    const impl_type* __impl;
  };
  
  
  template <typename Tp, typename Impl>
  inline bool
  operator==(const __snappy_vector_iterator<Tp,Impl>& x,
	     const __snappy_vector_iterator<Tp,Impl>& y)
  {
    return x.__impl == y.__impl && x.__pos == y.__pos;
  }
  
  template <typename Tp, typename Impl>
  inline bool
  operator!=(const __snappy_vector_iterator<Tp,Impl>& x,
	     const __snappy_vector_iterator<Tp,Impl>& y)
  {
    return x.__impl != y.__impl || x.__pos != y.__pos;
  }
  
  
  template <typename Tp, typename Impl>
  inline bool
  operator<(const __snappy_vector_iterator<Tp,Impl>& x,
	    const __snappy_vector_iterator<Tp,Impl>& y)
  {
    return ((x.__impl == y.__impl && x.__pos < y.__pos) || x.__impl < y.__impl);
  }
  
  template <typename Tp, typename Impl>
  inline bool
  operator>(const __snappy_vector_iterator<Tp,Impl>& x,
	    const __snappy_vector_iterator<Tp,Impl>& y)
  {
    return y < x;
  }
  
  
  template <typename Tp, typename Impl>
  inline bool
  operator<=(const __snappy_vector_iterator<Tp,Impl>& x,
	     const __snappy_vector_iterator<Tp,Impl>& y)
  {
    return ! (y < x);
  }
  
  template <typename Tp, typename Impl>
  inline bool
  operator>=(const __snappy_vector_iterator<Tp,Impl>& x,
	     const __snappy_vector_iterator<Tp,Impl>& y)
  {
    return ! (x < y);
  }

  template <typename Tp, typename Impl>
  inline ptrdiff_t
  operator-(const __snappy_vector_iterator<Tp,Impl>& x,
	    const __snappy_vector_iterator<Tp,Impl>& y)
  {
    return (x.__impl == y.__impl ? x.__pos - y.__pos : x.__impl - y.__impl);
  }
  
  template <typename Tp, typename Impl>
  inline __snappy_vector_iterator<Tp,Impl>
  operator+(ptrdiff_t __n, const __snappy_vector_iterator<Tp,Impl>& __x)
  {
    return __x + __n;
  }
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class snappy_vector_mapped : public __snappy_vector_base<Tp>
  {
  private:
    typedef __snappy_vector_base<Tp>        base_type;
    typedef snappy_vector_mapped<Tp, Alloc> self_type;
    
  public:
    typedef typename base_type::byte_type       byte_type;
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::value_type      value_type;
    
    typedef boost::filesystem::path            path_type;

    typedef __snappy_vector_iterator<Tp, self_type> const_iterator;
    typedef __snappy_vector_iterator<Tp, self_type>       iterator;
    
  private:
    typedef typename Alloc::template rebind<byte_type>::other file_allocator_type;
    typedef utils::snappy_file<file_allocator_type> file_type;

  public:
    snappy_vector_mapped() : __file(), __mutex_shift(0) {}
    snappy_vector_mapped(const path_type& path) : __mutex_shift(0) { open(path); }
    
  public:
    const_iterator begin() const { return const_iterator(0, *this); }
    const_iterator end() const { return const_iterator(size(), *this); }

    const value_type& operator[](size_type pos) const
    {
      const size_type buffer_pos = pos >> base_type::__buffer_shift;
      const size_type buffer_iter = pos - (buffer_pos * base_type::__buffer_size);
      
      const size_type cache_pos = buffer_pos & (__cache.size() - 1);
      
      lock_type lock(const_cast<spinlock_type&>(__mutex[cache_pos >> __mutex_shift].mutex));
      
      cache_type& cache = const_cast<cache_type&>(__cache[cache_pos]);
      if (cache.pos != buffer_pos) {
	cache.pos = buffer_pos;
	
	__file.read((byte_type*) &(*cache.buffer.begin()), base_type::__buffer_size * sizeof(value_type), buffer_pos * base_type::__buffer_size * sizeof(value_type));
      }
      
      return cache.buffer[buffer_iter];
    }
    
  public:
    void open(const path_type& path)
    {
      __file.open(path);
      
      // determine the cache size
      const size_type cache_byte_size = utils::bithack::max(utils::bithack::next_largest_power2(__file.size_compressed() >> 6),
							    size_type(1024 * 64));
      const size_type cache_size_raw = cache_byte_size / (base_type::__buffer_size * sizeof(value_type));
      const size_type cache_size_power2 = utils::bithack::branch(utils::bithack::is_power2(cache_size_raw),
								 cache_size_raw,
								 utils::bithack::next_largest_power2(cache_size_raw));
      
      __cache.clear();
      __cache.reserve(cache_size_power2);
      __cache.resize(cache_size_power2);
      
      __mutex_shift = utils::bithack::bit_count(cache_size_power2 - 1) - 3;
    }

    void write(const path_type& path) const
    {
      __file.write(path);
    }
    
    void close() { clear(); }
    void clear()
    {  
      __file.clear();
      __cache.clear();
      __mutex_shift = 0;
    }
    
    void swap(snappy_vector_mapped& x)
    {
      __file.swap(x.__file);
      __cache.swap(x.__cache);
      std::swap(__mutex_shift, x.__mutex_shift);
    }
    
    static bool exists(const path_type& path)
    {
      return file_type::exists(path);
    }
    
    size_type size() const { return __file.size() / sizeof(value_type); }
    bool empty() const { return __file.empty(); }
    path_type path() const { return __file.path(); }
    bool is_open() const { return __file.is_open(); }
    
    uint64_t size_bytes() const { return __file.size_bytes(); }
    uint64_t size_compressed() const { return __file.size_compressed(); }
    uint64_t size_cache() const { return __file.size_cache(); + sizeof(cache_type) * __cache.size(); }
    
  private:
    typedef utils::spinlock             spinlock_type;
    typedef spinlock_type::scoped_lock  lock_type;
    
    struct mutex_type
    {
      mutex_type() : mutex() {}
      mutex_type(const mutex_type& x) : mutex(){}
      mutex_type& operator=(const mutex_type& x) { return *this; }
      spinlock_type mutex;
    };
    
    typedef boost::array<mutex_type, 8> mutex_set_type;
    typedef boost::array<value_type, base_type::__buffer_size> buffer_type;

    struct cache_type
    {
      cache_type() : pos(size_type(-1)), buffer() {}

      size_type   pos;
      buffer_type buffer;
    };
    typedef typename Alloc::template rebind<cache_type>::other  cache_allocator_type;
    typedef std::vector<cache_type, cache_allocator_type> cache_set_type;
    
  private:
    file_type __file;

    cache_set_type __cache;
    mutex_set_type __mutex;
    size_type      __mutex_shift;
  };
  
};

namespace std
{
  template <typename T, typename A>
  inline
  void swap(utils::snappy_vector_mapped<T,A>& x, utils::snappy_vector_mapped<T,A>& y)
  {
    return x.swap(y);
  }
};

#endif
