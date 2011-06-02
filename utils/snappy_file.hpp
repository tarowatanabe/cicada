// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SNAPPY_FILE__HPP__
#define __UTILS__SNAPPY_FILE__HPP__ 1

#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>

#include <utils/vertical_coded_vector.hpp>
#include <utils/map_file.hpp>
#include <utils/repository.hpp>
#include <utils/filesystem.hpp>
#include <utils/bithack.hpp>
#include <utils/spinlock.hpp>

#include <snappy.h>

namespace utils
{
  template <typename Alloc=std::allocator<char> >
  class snappy_file
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef char      char_type;
    typedef char      byte_type;
    typedef uint64_t  off_type;

    typedef boost::filesystem::path path_type;
    
  private:
    typedef typename Alloc::template rebind<off_type>::other  off_vector_allocator_type;
    typedef typename Alloc::template rebind<byte_type>::other data_vector_allocator_type;
    
    typedef utils::vertical_coded_vector_mapped<off_type, off_vector_allocator_type> off_vector_type;
    typedef utils::map_file<byte_type, data_vector_allocator_type>                   data_vector_type;

    static const size_type __buffer_size = 4096;
    static const size_type __buffer_mask = 4096 - 1;
    static const size_type __buffer_shift = utils::bithack::static_bit_count<__buffer_mask>::value;

  public:
    snappy_file() : __size(0), __data(), __offset(), __mutex_shift(0) {}
    snappy_file(const path_type& path) : __size(0), __data(), __offset(), __mutex_shift(0) { open(path); }

  public:
    size_type read(byte_type* dest, size_type size, off_type offset) const
    {
      if (offset >= __size) return 0;
      
      for (size_type iter = 0; iter != size; /**/) {
	const size_type pos = offset + iter;
	const size_type buffer_pos = pos >> __buffer_shift;
	
	const size_type read_offset = pos - (buffer_pos * __buffer_size);
	const size_type read_size = utils::bithack::min(__buffer_size - read_offset, size - iter);
	
	const size_type cache_pos = buffer_pos & (__cache.size() - 1);
	
	lock_type lock(const_cast<spinlock_type&>(__mutex[cache_pos >> __mutex_shift].mutex));
	
	cache_type& cache = const_cast<cache_type&>(__cache[cache_pos]);
	if (cache.pos != buffer_pos) {
	  cache.pos = buffer_pos;
	  
	  const off_type first = __offset[buffer_pos];
	  const off_type last  = __offset[buffer_pos + 1];
	  
	  snappy::RawUncompress(__data.begin() + first, last - first, &(*cache.buffer.begin()));
	}
	
	std::copy(cache.buffer.begin() + read_offset, cache.buffer.begin() + read_offset + read_size, dest);
	
	iter += read_size;
	dest += read_size;
      }
      
      return size;
    }

  public:
    void open(const path_type& path)
    {
      typedef utils::repository repository_type;

      repository_type rep(path, repository_type::read);
      
      __data.open(rep.path("data"));
      __offset.open(rep.path("offsets"));
      
      repository_type::const_iterator iter = rep.find("size");
      if (iter == rep.end())
	throw std::runtime_error("no size?");
      __size = boost::lexical_cast<size_type>(iter->second);
      
      // determine the cache size...
      const size_type cache_byte_size = utils::bithack::max(size_type(utils::bithack::next_largest_power2(__data.size() >> 7)),
							    size_type(1024 * 64));
      const size_type cache_size_raw = cache_byte_size / __buffer_size;
      const size_type cache_size_power2 = utils::bithack::branch(utils::bithack::is_power2(cache_size_raw),
								 cache_size_raw,
								 size_type(utils::bithack::next_largest_power2(cache_size_raw)));
      
      __cache.clear();
      __cache.reserve(cache_size_power2);
      __cache.resize(cache_size_power2);
      
      __mutex_shift = utils::bithack::bit_count(cache_size_power2 - 1) - 3;
    }
    
    void write(const path_type& __path) const
    {
      if (__path == path()) return;
      
      typedef utils::repository repository_type;
      
      repository_type rep(__path, repository_type::write);
      
      __data.write(rep.path("data"));
      __offset.write(rep.path("offsets"));
      
      std::ostringstream stream_size;
      stream_size << __size;
      
      rep["type"] = "snappy";
      rep["size"] = stream_size.str();
      rep["block"] = "4096";
    }

    
    void close() { clear(); }
    void clear()
    {  
      __size = 0;
      __data.clear();
      __offset.clear();
      __cache.clear();
      __mutex_shift = 0;
    }

    void swap(snappy_file& x)
    {
      std::swap(__size, x.__size);
      __data.swap(x.__data);
      __offset.swap(x.__offset);
      __cache.swap(x.__cache);
      std::swap(__mutex_shift, x.__mutex_shift);
    }

    static bool exists(const path_type& path)
    {
      if (! utils::repository::exists(path)) return false;
      if (! off_vector_type::exists(path / "offsets")) return false;
      if (! data_vector_type::exists(path / "data")) return false;
      
      return true;
    }
    
    size_type size() const { return __size; }
    bool empty() const { return __size == 0; }
    path_type path() const { return __data.path().parent_path(); }
    bool is_open() const { return __data.is_open(); }
    
    uint64_t size_bytes() const { return __size; }
    uint64_t size_compressed() const { return __data.size_compressed() + __offset.size_compressed(); }
    uint64_t size_cache() const { return __data.size_cache() + __offset.size_cache() + sizeof(cache_type) * __cache.size(); }
    
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
    typedef boost::array<byte_type, 4096> buffer_type;
    
    struct cache_type
    {
      cache_type() : pos(size_type(-1)), buffer() {}
      
      size_type   pos;
      buffer_type buffer;
    };
    typedef typename Alloc::template rebind<cache_type>::other  cache_allocator_type;
    typedef std::vector<cache_type, cache_allocator_type> cache_set_type;
    
  private:
    size_type        __size;
    data_vector_type __data;
    off_vector_type  __offset;
    
    cache_set_type   __cache;
    mutex_set_type   __mutex;
    size_type        __mutex_shift;
  };
  
};

namespace std
{
  template <typename A>
  inline
  void swap(utils::snappy_file<A>& x, utils::snappy_file<A>& y)
  {
    return x.swap(y);
  }
};

#endif
