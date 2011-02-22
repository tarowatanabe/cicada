// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__VERTICAL_CODED_VECTOR__HPP__
#define __UTILS__VERTICAL_CODED_VECTOR__HPP__ 1

#include <stdint.h>

#include <vector>
#include <stdexcept>
#include <sstream>
#include <iterator>

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/thread.hpp>

#include <utils/atomicop.hpp>
#include <utils/bithack.hpp>
#include <utils/repository.hpp>
#include <utils/map_file.hpp>
#include <utils/array_power2.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/filesystem.hpp>


namespace utils
{
  
  template <typename Tp, typename Impl>
  struct __vertical_coded_vector_iterator
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Tp        value_type;
    typedef Impl      impl_type;
    
    typedef __vertical_coded_vector_iterator<Tp,Impl> self_type;

    typedef std::random_access_iterator_tag   iterator_category;
    
    __vertical_coded_vector_iterator(size_type pos, const impl_type* impl)
      : __pos(pos), __impl(impl) {}
    __vertical_coded_vector_iterator()
      : __pos(), __impl() {}
    
    value_type operator*() const { return __impl->operator[](__pos); }
    
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
  operator==(const __vertical_coded_vector_iterator<Tp,Impl>& x,
	     const __vertical_coded_vector_iterator<Tp,Impl>& y)
  {
    return x.__impl == y.__impl && x.__pos == y.__pos;
  }
  
  template <typename Tp, typename Impl>
  inline bool
  operator!=(const __vertical_coded_vector_iterator<Tp,Impl>& x,
	     const __vertical_coded_vector_iterator<Tp,Impl>& y)
  {
    return x.__impl != y.__impl || x.__pos != y.__pos;
  }
  
  
  template <typename Tp, typename Impl>
  inline bool
  operator<(const __vertical_coded_vector_iterator<Tp,Impl>& x,
	    const __vertical_coded_vector_iterator<Tp,Impl>& y)
  {
    return ((x.__impl == y.__impl && x.__pos < y.__pos) || x.__impl < y.__impl);
  }
  
  template <typename Tp, typename Impl>
  inline bool
  operator>(const __vertical_coded_vector_iterator<Tp,Impl>& x,
	    const __vertical_coded_vector_iterator<Tp,Impl>& y)
  {
    return y < x;
  }
  
  
  template <typename Tp, typename Impl>
  inline bool
  operator<=(const __vertical_coded_vector_iterator<Tp,Impl>& x,
	     const __vertical_coded_vector_iterator<Tp,Impl>& y)
  {
    return ! (y < x);
  }
  
  template <typename Tp, typename Impl>
  inline bool
  operator>=(const __vertical_coded_vector_iterator<Tp,Impl>& x,
	     const __vertical_coded_vector_iterator<Tp,Impl>& y)
  {
    return ! (x < y);
  }

  template <typename Tp, typename Impl>
  inline ptrdiff_t
  operator-(const __vertical_coded_vector_iterator<Tp,Impl>& x,
	    const __vertical_coded_vector_iterator<Tp,Impl>& y)
  {
    return (x.__impl == y.__impl ? x.__pos - y.__pos : x.__impl - y.__impl);
  }

  template <typename Tp, typename Impl>
  inline __vertical_coded_vector_iterator<Tp,Impl>
  operator+(ptrdiff_t __n, const __vertical_coded_vector_iterator<Tp,Impl>& __x)
  {
    return __x + __n;
  }
  
  struct __vertical_coded_vector_base
  {
    typedef uint8_t   byte_type;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef uint64_t  off_type;
    
    typedef boost::filesystem::path path_type;
    
    template <typename Compressed, typename Off>
    off_type read(const Compressed& compressed, const Off& off, const size_type pos, const size_type iter) const
    {
      const size_type pos_first = off[iter];
      const size_type pos_last  = off[iter + 1];
      const size_type pos_access = pos_first + pos;
      
      return (pos_access >= pos_last
	      ? off_type(-1)
	      : (iter == off.size() - 2
		 ? off_type(compressed[pos_access])
		 : (last(compressed, off, pos, iter + 1) << 8) | off_type(compressed[pos_access])));
    }
    
    template <typename Compressed, typename Off>
    off_type last(const Compressed& compressed, const Off& off, const size_type pos, const size_type iter) const
    {
      const size_type pos_first = off[iter];
      const size_type pos_last  = off[iter + 1];
      
      if (iter == off.size() - 2)
	return (pos > 0xff ? off_type(pos_last - pos_first) : upper_bound(compressed, pos_first, pos_last, pos) - pos_first) - 1;
      else {
	const off_type pos_bound_first = read(compressed, off, pos >> 8,       iter + 1);
	const off_type pos_bound_last  = read(compressed, off, (pos >> 8) + 1, iter + 1);

	const size_type pos_bound_first_mask = size_type(pos_bound_first == off_type(-1)) - 1;
	const size_type pos_bound_last_mask  = size_type(pos_bound_last == off_type(-1)) - 1;

	const size_type pos_access_first = ((~pos_bound_first_mask) & pos_last) | (pos_bound_first_mask & (pos_first + pos_bound_first));
	const size_type pos_access_last  = ((~pos_bound_last_mask)  & pos_last) | (pos_bound_last_mask  & (pos_first + pos_bound_last));
	
	//const size_type pos_access_first = (pos_bound_first == off_type(-1) ? pos_last : pos_first + pos_bound_first);
	//const size_type pos_access_last  = (pos_bound_last  == off_type(-1) ? pos_last : pos_first + pos_bound_last);
	
	return upper_bound(compressed, pos_access_first, pos_access_last, pos & 0xff) - pos_first - 1;
      }
    }
    
    template <typename Compressed>
    off_type upper_bound(const Compressed& compressed, size_type first, size_type last, const byte_type& value) const
    {
      size_type length = last - first;
      if (length < 128) {
	typename Compressed::const_iterator iter = compressed.begin() + first;
	for (/**/; first != last && (*iter) <= value; ++ first, ++ iter);
	return first;
      } else {
	typename Compressed::const_iterator __first = compressed.begin() + first;
	typename Compressed::const_iterator __middle;
	
	while (length > 0) {
	  const size_type half = length >> 1;
	  const size_type middle = first + half;
	  __middle = __first + half;
	  
	  if (value < *__middle)
	    length = half;
	  else {
	    first = middle + 1;
	    __first = __middle + 1;
	    length = length - half - 1;
	  }
	}
	return first;
      }
    }
  };


  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class vertical_coded_vector_mapped : public __vertical_coded_vector_base
  {
  public:
    typedef Tp        value_type;
    
    typedef __vertical_coded_vector_base base_type;
    
    typedef base_type::byte_type       byte_type;
    typedef base_type::size_type       size_type;
    typedef base_type::difference_type difference_type;
    typedef base_type::off_type        off_type;
    typedef base_type::path_type       path_type;
    
  private:
    typedef typename Alloc::template rebind<byte_type>::other  byte_allocator_type;
    typedef typename Alloc::template rebind<off_type>::other   off_allocator_type;
    
    typedef utils::map_file<byte_type, byte_allocator_type>   compressed_vector_type;
    typedef utils::map_file<off_type, off_allocator_type>     off_vector_type;

  public:
    typedef vertical_coded_vector_mapped<Tp, Alloc>         self_type;
    typedef __vertical_coded_vector_iterator<Tp, self_type> const_iterator;
    typedef __vertical_coded_vector_iterator<Tp, self_type> iterator;

  private:
    struct Cache
    {
     typedef int64_t value_type;
      
      Cache() : value(value_type(-1)) {}
      
      volatile value_type value;
    };
    typedef Cache cache_type;
    typedef typename Alloc::template rebind<cache_type>::other   cache_allocator_type;
    typedef std::vector<cache_type, cache_allocator_type> cache_set_type;
    
  public:
    vertical_coded_vector_mapped(const path_type& path) { open(path); }
    vertical_coded_vector_mapped() {}
    
  public:
    
    bool empty() const { return compressed.empty(); }
    size_type size() const { return empty() ? size_type(0) : size_type(off[1]); }
    path_type path() const { return compressed.path().parent_path(); }
    bool is_open() const { return compressed.is_open(); }
    
    uint64_t size_bytes() const { return size() * sizeof(Tp); }
    uint64_t size_compressed() const { return compressed.size_compressed() + off.size_compressed(); }
    uint64_t size_cache() const { return __cache.size() * sizeof(cache_type); }
    
    void close() { clear(); }
    void clear()
    {
      compressed.clear();
      off.clear();
      __cache.clear();
    }
    
    const_iterator begin() const { return const_iterator(size_type(0), this); }
    const_iterator end() const { return const_iterator(size(), this); }

    value_type front() const 
    {
      return operator[](0);
    }
    
    value_type back() const 
    {
      return operator[](size() - 1);
    }
    
    value_type operator[](size_type pos) const
    {
      cache_set_type& caches = const_cast<cache_set_type&>(__cache);

      const uint64_t mask_cache = caches.size() - 1;
      
      cache_type cache;
      cache_type cache_new;
      
      cache.value = utils::atomicop::fetch_and_add(caches[pos & mask_cache].value, int64_t(0));
      
      uint64_t __pos   = (cache.value & __mask_pos);
      uint64_t __value = (cache.value & __mask_value);
      
      if (__pos == ((uint64_t(pos) << 32) & __mask_pos))
	return __value;
      
      __value = read(compressed, off, pos, 0);

      cache_new.value = ((uint64_t(pos) << 32) & __mask_pos) | (uint64_t(__value) & __mask_value);
      
      utils::atomicop::compare_and_swap(caches[pos & mask_cache].value, cache.value, cache_new.value);
      
      return __value;
    }
    
  public:
    static bool exists(const path_type& path)
    {
      if (! utils::repository::exists(path)) return false;
      if (! compressed_vector_type::exists(path / "data")) return false;
      if (! off_vector_type::exists(path / "offsets")) return false;
      
      return true;
    }

    void open(const path_type& path)
    {
      typedef utils::repository repository_type;
      
      clear();

      repository_type rep(path, repository_type::read);
      compressed.open(rep.path("data"));
      off.open(rep.path("offsets"));
      
      repository_type::const_iterator titer = rep.find("type");
      if (titer == rep.end())
	throw std::runtime_error("no type...");
      if (titer->second != "vertical-coded")
	throw std::runtime_error("not a vertica-coded");
      
      const size_type cache_size = std::max(size_type(utils::bithack::next_largest_power2(size() >> 7)),
					    size_type(1024 * 32));
      
      __cache.reserve(cache_size);
      __cache.resize(cache_size, cache_type());
      
      __mask_pos = (~uint64_t(__cache.size() - 1)) << 32;
      __mask_value = ~__mask_pos;
    }

    void write(const path_type& file) const
    {
      if (path() == file) return;
      
      // remove first...
      if (boost::filesystem::exists(file) && ! boost::filesystem::is_directory(file))
	boost::filesystem::remove_all(file);
      
      // create directory
      if (! boost::filesystem::exists(file))
	boost::filesystem::create_directories(file);
      
      // remove all the files...
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(file); iter != iter_end; ++ iter)
	boost::filesystem::remove_all(*iter);
      
      // copy all...
      for (boost::filesystem::directory_iterator iter(path()); iter != iter_end; ++ iter)
	utils::filesystem::copy_files(*iter, file);
    }
    
  private:
    compressed_vector_type compressed;
    off_vector_type        off;
    
    cache_set_type __cache;
    uint64_t       __mask_pos;
    uint64_t       __mask_value;
  };
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class vertical_coded_vector : public __vertical_coded_vector_base
  {
  public:
    typedef Tp        value_type;
    
    typedef __vertical_coded_vector_base base_type;
    
    typedef base_type::byte_type       byte_type;
    typedef base_type::size_type       size_type;
    typedef base_type::difference_type difference_type;
    typedef base_type::off_type        off_type;
    typedef base_type::path_type       path_type;
    
  private:
    typedef typename Alloc::template rebind<byte_type>::other  byte_allocator_type;
    typedef typename Alloc::template rebind<value_type>::other value_allocator_type;
    typedef typename Alloc::template rebind<off_type>::other   off_allocator_type;
    
    typedef std::vector<byte_type, byte_allocator_type>   compressed_vector_type;
    typedef std::vector<off_type, off_allocator_type>     off_vector_type;
    typedef std::vector<value_type, value_allocator_type> raw_vector_type;

  public:
    typedef vertical_coded_vector<Tp, Alloc>         self_type;
    typedef __vertical_coded_vector_iterator<Tp, self_type> const_iterator;
    typedef __vertical_coded_vector_iterator<Tp, self_type> iterator;

  public:
    vertical_coded_vector() {}
    
  public:

    bool empty() const { return compressed.empty() && raw.empty(); }
    size_type size() const { return (off.empty() ? raw.size() : size_type(off[1])); }
    
    uint64_t size_bytes() const { return size() * sizeof(Tp); }
    uint64_t size_compressed() const { return compressed.size() * sizeof(byte_type) + off.size() * sizeof(off_type); }
    uint64_t size_cache() const { return 0; }
        
    void clear()
    {
      compressed.clear();
      off.clear();
      raw.clear();
    }

    const_iterator begin() const { return const_iterator(size_type(0), this); }
    const_iterator end() const { return const_iterator(size(), this); }

    void push_back(const value_type& x)
    {
      if (! compressed.empty())
	throw std::runtime_error("modification not supported for compressed vector...");
      raw.push_back(x);
    }
    
    template <typename Iterator, typename InsertIterator>
    void insert(Iterator pos, InsertIterator first, InsertIterator last)
    {
      if (! compressed.empty())
	throw std::runtime_error("modification not supported for compressed vector...");
      
      raw.insert(pos, first, last);
    }
    
    void resize(const size_type new_size)
    {
      if (! compressed.empty())
	throw std::runtime_error("modification not supported for compressed vector...");
      raw.resize(new_size);
    }
    
    void resize(const size_type new_size, const value_type& x)
    {
      if (! compressed.empty())
	throw std::runtime_error("modification not supported for compressed vector...");
      
      raw.resize(new_size, x);
    }

    value_type& front() 
    {
      if (! compressed.empty())
	throw std::runtime_error("modification not supported for compressed vector...");
      return raw.front();
    }

    value_type& back() 
    {
      if (! compressed.empty())
	throw std::runtime_error("modification not supported for compressed vector...");
      return raw.back();
    }
    
    value_type front() const 
    {
      return (compressed.empty() ? raw.front() : operator[](0));
    }
    
    value_type back() const 
    {
      return (compressed.empty() ? raw.back() : operator[](size() - 1));
    }
    
    value_type& operator[](const size_type pos)
    {
      if (! compressed.empty())
	throw std::runtime_error("modification not supported for compressed vector...");
      return raw[pos];
    }
    value_type operator[](size_type pos) const
    {
      return (compressed.empty() ? raw[pos] : read(compressed, off, pos, 0));
    }
    
  public:
    
    void build()
    {
      typedef raw_vector_type inverted_vector_type;
      
      compressed.clear();
      off.clear();
      
      if (raw.empty()) return;
      
      inverted_vector_type inverted;
      inverted_vector_type inverted_new;
      
      bool finished = false;
      
      for (size_type i = 0; i < raw.size(); ++ i) {
	const byte_type  value  = raw[i] & 0xff;
	const value_type invert = raw[i] >> 8;
	
	compressed.push_back(value);
	
	if (inverted.empty() || value_type(inverted.size() - 1) != invert)
	  inverted.resize(invert + 1, i);
      }
      off.push_back(0);
      off.push_back(compressed.size());
      
      if (inverted.back() <= 0xff)
	finished = true;
      
      while (! finished) {
	
	for (size_type i = 0; i < inverted.size(); ++ i) {
	  const byte_type  value  = inverted[i] & 0xff;
	  const value_type invert = inverted[i] >> 8;
	  
	  compressed.push_back(value);
	  
	  if (inverted_new.empty() || value_type(inverted_new.size() - 1) != invert)
	    inverted_new.resize(invert + 1, i);
	}
	
	inverted_new.swap(inverted);
	inverted_new.clear();
	
	off.push_back(compressed.size());
	
	if (inverted.back() <= 0xff)
	  finished = true;
      }
      
      compressed.insert(compressed.end(), inverted.begin(), inverted.end());
      off.push_back(compressed.size());
    }
    
  public:
    void write(const path_type& path) const
    {
      typedef utils::repository repository_type;

      if (compressed.empty())
	const_cast<vertical_coded_vector&>(*this).build();
      
      repository_type rep(path, repository_type::write);
      dump_file(rep.path("data"), compressed);
      dump_file(rep.path("offsets"), off);
      
      std::ostringstream stream_integral_size;
      std::ostringstream stream_size;
      std::ostringstream stream_coded_size;
      
      stream_integral_size << sizeof(value_type);
      stream_size << size();
      stream_coded_size << off.size() - 1;
      
      rep["size"] = stream_size.str();
      rep["integral-size"] = stream_integral_size.str();
      rep["code-size"] = stream_coded_size.str();
      rep["type"] = "vertical-coded";
    }
    
  private:
    template <typename _Path, typename _Data>
    inline
    void dump_file(const _Path& file, const _Data& data) const
    {
      boost::iostreams::filtering_ostream os;
#if BOOST_FILESYSTEM_VERSION == 2
      os.push(boost::iostreams::file_sink(file.native_file_string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
#else
      os.push(boost::iostreams::file_sink(file.string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
#endif
      
      const int64_t file_size = sizeof(typename _Data::value_type) * data.size();
      for (int64_t offset = 0; offset < file_size; offset += 1024 * 1024)
	if (! os.write(((char*) &(*data.begin())) + offset, std::min(int64_t(1024 * 1024), file_size - offset)))
	  throw std::runtime_error("vertical_coded_vector::write()");
    }
    
  private:
    compressed_vector_type compressed;
    off_vector_type        off;
    raw_vector_type        raw;
  };
  
};

#endif
