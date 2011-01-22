// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__PACKED_VECTOR__HPP__
#define __UTILS__PACKED_VECTOR__HPP__ 1

// byte packed vector with random access support


#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>

#include <boost/filesystem.hpp>

#include <utils/succinct_vector.hpp>
#include <utils/map_file.hpp>
#include <utils/repository.hpp>
#include <utils/filesystem.hpp>

namespace utils
{

  template <typename Tp, size_t Size>
  struct __packed_vector_base {};

  // use of 4-bit alignment...?
  template <typename Tp>
  struct __packed_vector_base<Tp, 1>
  {
    typedef uint8_t   byte_type;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef Tp        value_type;
    typedef uint8_t   uvalue_type;
    
    size_type byte_size(const value_type& x) const
    {
      return 1 + bool(uvalue_type(x) & 0xf0);
    }
    
    template <typename Data, typename Index>
    size_type encode(Data& data, size_type& coded_size, Index& index, const value_type& value)
    {
      const size_type value_size = byte_size(value);
      
      if (coded_size & 0x01) {
	switch (value_size) {
	case 2:
	  data.back() |= (value & 0x0f) << 4;
	  data.push_back((value >> 4) & 0x0f);
	  break;
	case 1:
	  data.back() |= (value & 0x0f) << 4;
	  break;
	}
      } else
	data.push_back(value);
      
      index.set(index.size(), value_size == 2);
      coded_size += value_size;
      
      return coded_size;
    }
    
    template <typename Data, typename Index>
    value_type decode(const Data& data, const Index& index, size_type pos) const
    {
      const size_type pos_first = (pos == 0 ? size_type(0) : pos - 1 + index.rank(pos - 1, true) + 1);
      const size_type pos_last  = pos + index.rank(pos, true) + 1;
      
      const size_type pos_byte = pos_first >> 1;
      if (pos_first & 0x01) {
	switch (pos_last - pos_first) {
	case 2: return ((data[pos_byte] >> 4) & 0x0f) | ((data[pos_byte + 1] << 4) & 0xf0);
	case 1: return (data[pos_byte] >> 4) & 0x0f;
	}
      } else {
	switch (pos_last - pos_first) {
	case 2: return data[pos_byte];
	case 1: return data[pos_byte] & 0x0f;
	}
      }
      
      // we will never reach here...
      return 0;
    }
  };
  
  template <typename Tp>
  struct __packed_vector_base<Tp, 2>
  {
    typedef uint8_t   byte_type;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;    
    typedef Tp        value_type;
    typedef uint16_t  uvalue_type;
    
    size_type byte_size(const value_type& x) const
    {
      return (1 + bool(uvalue_type(x) & 0xff00));
    }
    
    template <typename Data, typename Index>
    size_type encode(Data& data, size_type& coded_size, Index& index, const value_type& value)
    {
      const size_type value_size = byte_size(value);
      const size_type pos_first = data.size();
      
      data.resize(data.size() + value_size);
      
      switch (value_size) {
      case 2: data[pos_first + 1] = (uvalue_type(value) >> 8);
      case 1: data[pos_first + 0] = (uvalue_type(value));
      }
      
      index.set(index.size(), value_size == 2);
      coded_size += value_size;
      
      return coded_size;
    }
    
    template <typename Data, typename Index>
    value_type decode(const Data& data, const Index& index, size_type pos) const
    {
      const size_type pos_first = (pos == 0 ? size_type(0) : pos - 1 + index.rank(pos - 1, true) + 1);
      const size_type pos_last  = pos + index.rank(pos, true) + 1;
      
      const uvalue_type mask = 0xff;
      
      uvalue_type value = 0;
      switch (pos_last - pos_first) {
      case 2: value |= ((uvalue_type(data[pos_first + 1]) & mask) << 8);
      case 1: value |= ((uvalue_type(data[pos_first + 0]) & mask));
      }
      
      return value_type(value);
    }
  };
  
  template <typename Tp>
  struct __packed_vector_base<Tp, 4>
  {
    typedef uint8_t   byte_type;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;    
    typedef Tp        value_type;
    typedef uint32_t  uvalue_type;
    
    size_type byte_size(const value_type& x) const
    {
      return __byte_size((uvalue_type) x);
    }
    
    size_type __byte_size(uvalue_type x) const
    {
      return (1 
	      + bool(x & 0xffffff00)
	      + bool(x & 0xffff0000)
	      + bool(x & 0xff000000));
    }
    
    template <typename Data, typename Index>
    size_type encode(Data& data, size_type& coded_size, Index& index, const value_type& value)
    {
      const size_type value_size = byte_size(value);
      const size_type pos_first = data.size();
      
      data.resize(data.size() + value_size);
      
      switch (value_size) {
      case 4: data[pos_first + 3] = (uvalue_type(value) >> 24);
      case 3: data[pos_first + 2] = (uvalue_type(value) >> 16);
      case 2: data[pos_first + 1] = (uvalue_type(value) >> 8);
      case 1: data[pos_first + 0] = (uvalue_type(value));
      }
      
      index.set(index.size() + value_size - 1, true);
      coded_size += value_size;
      
      return coded_size;
    }
    
    template <typename Data, typename Index>
    value_type decode(const Data& data, const Index& index, size_type pos) const
    {
      const size_type pos_first = (pos == 0 ? size_type(0) : size_type(index.select(pos, true) + 1));
      const size_type pos_last = index.select(pos + 1, true) + 1;
      
      const uvalue_type mask = 0xff;
      
      uvalue_type value = 0;
      switch (pos_last - pos_first) {
      case 4: value |= ((uvalue_type(data[pos_first + 3]) & mask) << 24);
      case 3: value |= ((uvalue_type(data[pos_first + 2]) & mask) << 16);
      case 2: value |= ((uvalue_type(data[pos_first + 1]) & mask) << 8);
      case 1: value |= ((uvalue_type(data[pos_first + 0]) & mask));
      }
      
      return value_type(value);
    }
  };
  
  template <typename Tp>
  struct __packed_vector_base<Tp, 8>
  {
    typedef uint8_t   byte_type;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef Tp        value_type;
    typedef uint64_t  uvalue_type;
    
    size_type byte_size(const value_type& x) const
    {
      return __byte_size((uvalue_type) x);
    }
    
    size_type __byte_size(const value_type& x) const
    {
      return (1 
	      + bool(x & 0xffffffffffffff00ull)
	      + bool(x & 0xffffffffffff0000ull)
	      + bool(x & 0xffffffffff000000ull)
	      + bool(x & 0xffffffff00000000ull)
	      + bool(x & 0xffffff0000000000ull)
	      + bool(x & 0xffff000000000000ull)
	      + bool(x & 0xff00000000000000ull));
    }
    
    template <typename Data, typename Index>
    size_type encode(Data& data, size_type& coded_size, Index& index, const value_type& value)
    {
      const size_type value_size = byte_size(value);
      const size_type pos_first = data.size();
      
      data.resize(data.size() + value_size);
      
      switch (value_size) {
      case 8: data[pos_first + 7] = (uvalue_type(value) >> 56);
      case 7: data[pos_first + 6] = (uvalue_type(value) >> 48);
      case 6: data[pos_first + 5] = (uvalue_type(value) >> 40);
      case 5: data[pos_first + 4] = (uvalue_type(value) >> 32);
      case 4: data[pos_first + 3] = (uvalue_type(value) >> 24);
      case 3: data[pos_first + 2] = (uvalue_type(value) >> 16);
      case 2: data[pos_first + 1] = (uvalue_type(value) >> 8);
      case 1: data[pos_first + 0] = (uvalue_type(value));
      }
      
      index.set(index.size() + value_size - 1, true);
      coded_size += value_size;
      
      return coded_size;
    }
    
    template <typename Data, typename Index>
    value_type decode(const Data& data, const Index& index, size_type pos) const
    {
      const size_type pos_first = (pos == 0 ? size_type(0) : size_type(index.select(pos, true) + 1));
      const size_type pos_last = index.select(pos + 1, true) + 1;
      
      const uvalue_type mask = 0xff;
      
      uvalue_type value = 0;
      switch (pos_last - pos_first) {
      case 8: value |= ((uvalue_type(data[pos_first + 7]) & mask) << 56);
      case 7: value |= ((uvalue_type(data[pos_first + 6]) & mask) << 48);
      case 6: value |= ((uvalue_type(data[pos_first + 5]) & mask) << 40);
      case 5: value |= ((uvalue_type(data[pos_first + 4]) & mask) << 32);
      case 4: value |= ((uvalue_type(data[pos_first + 3]) & mask) << 24);
      case 3: value |= ((uvalue_type(data[pos_first + 2]) & mask) << 16);
      case 2: value |= ((uvalue_type(data[pos_first + 1]) & mask) << 8);
      case 1: value |= ((uvalue_type(data[pos_first + 0]) & mask));
      }
      
      return value_type(value);
    }
  };

  
  template <typename Tp, typename Impl>
  struct __packed_vector_iterator
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Tp        value_type;
    typedef Impl      impl_type;
    
    typedef __packed_vector_iterator<Tp,Impl> self_type;

    typedef Tp* pointer;
    typedef const Tp& reference;
    typedef std::random_access_iterator_tag   iterator_category;
    
    __packed_vector_iterator(size_type pos, const impl_type* impl)
      : __pos(pos), __impl(impl) { }
    __packed_vector_iterator()
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
  operator==(const __packed_vector_iterator<Tp,Impl>& x,
	     const __packed_vector_iterator<Tp,Impl>& y)
  {
    return x.__impl == y.__impl && x.__pos == y.__pos;
  }
  
  template <typename Tp, typename Impl>
  inline bool
  operator!=(const __packed_vector_iterator<Tp,Impl>& x,
	     const __packed_vector_iterator<Tp,Impl>& y)
  {
    return x.__impl != y.__impl || x.__pos != y.__pos;
  }
  
  
  template <typename Tp, typename Impl>
  inline bool
  operator<(const __packed_vector_iterator<Tp,Impl>& x,
	    const __packed_vector_iterator<Tp,Impl>& y)
  {
    return ((x.__impl == y.__impl && x.__pos < y.__pos) || x.__impl < y.__impl);
  }
  
  template <typename Tp, typename Impl>
  inline bool
  operator>(const __packed_vector_iterator<Tp,Impl>& x,
	    const __packed_vector_iterator<Tp,Impl>& y)
  {
    return y < x;
  }
  
  
  template <typename Tp, typename Impl>
  inline bool
  operator<=(const __packed_vector_iterator<Tp,Impl>& x,
	     const __packed_vector_iterator<Tp,Impl>& y)
  {
    return ! (y < x);
  }
  
  template <typename Tp, typename Impl>
  inline bool
  operator>=(const __packed_vector_iterator<Tp,Impl>& x,
	     const __packed_vector_iterator<Tp,Impl>& y)
  {
    return ! (x < y);
  }

  template <typename Tp, typename Impl>
  inline ptrdiff_t
  operator-(const __packed_vector_iterator<Tp,Impl>& x,
	    const __packed_vector_iterator<Tp,Impl>& y)
  {
    return (x.__impl == y.__impl ? x.__pos - y.__pos : x.__impl - y.__impl);
  }
  
  template <typename Tp, typename Impl>
  inline __packed_vector_iterator<Tp,Impl>
  operator+(ptrdiff_t __n, const __packed_vector_iterator<Tp,Impl>& __x)
  {
    return __x + __n;
  }
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class packed_vector_mapped : public __packed_vector_base<Tp, sizeof(Tp)>
  {
  private:
    typedef __packed_vector_base<Tp, sizeof(Tp)> base_type;
    typedef packed_vector_mapped<Tp, Alloc> self_type;
    
  public:
    typedef typename base_type::byte_type       byte_type;
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::value_type      value_type;
    typedef typename base_type::uvalue_type     uvalue_type;

    typedef boost::filesystem::path            path_type;

    typedef __packed_vector_iterator<Tp, self_type> const_iterator;
    typedef __packed_vector_iterator<Tp, self_type>       iterator;
    
  private:
    typedef typename Alloc::template rebind<byte_type>::other index_vector_allocator_type;
    typedef typename Alloc::template rebind<byte_type>::other data_vector_allocator_type;
    
    typedef utils::succinct_vector_mapped<index_vector_allocator_type> index_vector_type;
    typedef utils::map_file<byte_type, data_vector_allocator_type>  data_vector_type;
    
  public:
    packed_vector_mapped() : __size(0), __data(), __index() {}
    packed_vector_mapped(const path_type& path) { open(path); }
    
  public:
    const_iterator begin() const { return const_iterator(0, this); }
    const_iterator end() const { return const_iterator(__size, this); }

    value_type operator[](size_type pos) const
    {
      return base_type::decode(__data, __index, pos);
    }
    
    size_type size() const { return __size; }
    bool empty() const { return __size == 0; }
    path_type path() const { return __index.path().parent_path(); }
    bool is_open() const { return __index.is_open(); }
    
    uint64_t size_bytes() const { return sizeof(Tp) * __size; }
    uint64_t size_compressed() const { return __data.size_compressed() + __index.size_compressed(); }
    uint64_t size_cache() const { return __data.size_cache() + __index.size_cache(); }

    void close() { clear(); }
    void clear()
    {
      __data.clear();
      __index.clear();
      __size = 0;
    }

    static bool exists(const path_type& path)
    {
      if (! utils::repository::exists(path)) return false;
      if (! index_vector_type::exists(path / "index")) return false;
      if (! data_vector_type::exists(path / "data")) return false;
      
      return true;
    }
    
    void open(const path_type& path)
    {
      typedef utils::repository repository_type;
      
      repository_type rep(path, repository_type::read);
      
      __data.open(rep.path("data"));
      __index.open(rep.path("index"));
      
      repository_type::const_iterator iter = rep.find("size");
      if (iter == rep.end())
	throw std::runtime_error("no size?");
      __size = atoll(iter->second.c_str());
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

  public:
    size_type         __size;
    data_vector_type  __data;
    index_vector_type __index;
  };
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class packed_vector : public __packed_vector_base<Tp, sizeof(Tp)>
  {
  private:
    typedef __packed_vector_base<Tp, sizeof(Tp)> base_type;
    typedef packed_vector<Tp, Alloc> self_type;
    
  public:
    typedef typename base_type::byte_type       byte_type;
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::value_type      value_type;
    typedef typename base_type::uvalue_type     uvalue_type;

    typedef boost::filesystem::path            path_type;

    typedef __packed_vector_iterator<Tp, self_type> const_iterator;
    typedef __packed_vector_iterator<Tp, self_type>       iterator;
    
  private:
    typedef typename Alloc::template rebind<byte_type>::other index_vector_allocator_type;
    typedef typename Alloc::template rebind<byte_type>::other data_vector_allocator_type;
    
    typedef utils::succinct_vector<index_vector_allocator_type> index_vector_type;
    typedef std::vector<byte_type, data_vector_allocator_type>  data_vector_type;

  public:
    packed_vector() : __size(0), __size_coded(0), __data(), __index() {}
    
  public:
    
    const_iterator begin() const { return const_iterator(0, this); }
    const_iterator end() const { return const_iterator(__size, this); }

    size_type push_back(const value_type& x)
    {
      ++ __size;
      return base_type::encode(__data, __size_coded, __index, x);
    }
    
    template <typename Iterator>
    size_type insert(Iterator first, Iterator last)
    {
      for (/**/; first != last; ++ first, ++ __size)
	base_type::encode(__data, __size_coded, __index, *first);
      return __size_coded;
    }
    
    value_type operator[](size_type pos) const
    {
      return base_type::decode(__data, __index, pos);
    }
    
    size_type size() const { return __size; }
    bool empty() const { return __size == 0; }

    uint64_t size_bytes() const { return sizeof(Tp) * __size; }
    uint64_t size_compressed() const { return __data.size() * sizeof(byte_type) + __index.size_compressed(); }
    uint64_t size_cache() const { return __index.size_cache(); }
    
    void clear()
    {
      __data.clear();
      __index.clear();
      __size = 0;
      __size_coded = 0;
    }
    
    void build()
    {
      __index.build();
    }

    void write(const path_type& path) const
    {
      typedef utils::repository repository_type;
      
      repository_type rep(path, repository_type::write);
      
      dump_file(rep.path("data"), __data);
      __index.write(rep.path("index"));
      
      std::ostringstream stream_size;
      std::ostringstream stream_integral_size;
      stream_size << __size;
      stream_integral_size << sizeof(Tp);
      
      rep["size"] = stream_size.str();
      rep["integral-size"] = stream_integral_size.str();
      rep["type"] = "packed";
    }
    
  private:
    template <typename _Path, typename _Data>
    inline
    void dump_file(const _Path& file, const _Data& data) const
    {
      boost::iostreams::filtering_ostream os;
      os.push(boost::iostreams::file_sink(file.native_file_string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
      
      const int64_t file_size = sizeof(typename _Data::value_type) * data.size();
      for (int64_t offset = 0; offset < file_size; offset += 1024 * 1024)
	if (! os.write(((char*) &(*data.begin())) + offset, std::min(int64_t(1024 * 1024), file_size - offset)))
	  throw std::runtime_error("packed_vector::write()");
    }
    
  public:
    size_type         __size;
    size_type         __size_coded;
    data_vector_type  __data;
    index_vector_type __index;
  };
  
};

#endif
