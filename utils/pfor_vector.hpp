// -*- mode: c++ -*-

#ifndef __UTILS__PFOR_VECTOR__HPP__
#define __UTILS__PFOR_VECTOR__HPP__ 1

#include <stdint.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>

#include <boost/filesystem.hpp>

#include <utils/repository.hpp>
#include <utils/map_file.hpp>
#include <utils/bithack.hpp>
#include <utils/bitpack.hpp>

#include <codec/pfor_codec.hpp>

namespace utils
{

  template <typename Tp>
  struct __pfor_value_proxy
  {
    typedef Tp value_type;
    typedef Tp encoded_value_type;
    typedef Tp decoded_value_type;
    
    static inline
    Tp encode(const Tp& x) { return x; }
    
    static inline
    Tp decode(const Tp& x) { return x; }
  };
  
  template <>
  struct __pfor_value_proxy<float>
  {
    typedef float    value_type;
    typedef uint32_t encoded_value_type;
    typedef float    decoded_value_type;

    static const encoded_value_type mask = (encoded_value_type(1) << 31);
    
    static inline
    value_type decode(const encoded_value_type& x)
    {
      return __cast(x & mask ? (x & ~mask) : ~x);
    }
    
    static inline
    encoded_value_type encode(const value_type& x)
    {
      const encoded_value_type uvalue = ((encoded_value_type&) x);
      return (uvalue & mask ? ~uvalue : (uvalue | mask));
    }
    
    static inline
    value_type __cast(const encoded_value_type& x)
    {
      return (value_type&) x;
    }
  };
  
  template <>
  struct __pfor_value_proxy<double>
  {
    typedef float    value_type;
    typedef uint64_t encoded_value_type;
    typedef float    decoded_value_type;

    static const encoded_value_type mask = (encoded_value_type(1) << 63);

    static inline
    value_type decode(const encoded_value_type& x)
    { 
      return __cast(x & mask ? (x & ~mask) : ~x);
    }
    
    static inline
    encoded_value_type encode(const value_type& x)
    {
      const encoded_value_type uvalue = ((encoded_value_type&) x);
      return (uvalue & mask ? ~uvalue : (uvalue | mask));
    }
    
    static inline
    value_type __cast(const encoded_value_type& x)
    {
      return (value_type&) x;
    }
  };

  template <typename Tp, typename Impl>
  struct __pfor_vector_iterator
  {
  private:
    typedef Impl impl_type;

  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef int64_t   off_type;
    typedef uint32_t  pos_type;
    typedef Tp        value_type;
    
    typedef std::random_access_iterator_tag   iterator_category;
    
  private:
    typedef __pfor_vector_iterator<Tp, Impl> self_type;
    
  public:
    
    __pfor_vector_iterator(size_type pos, const impl_type* impl)
      : __pos(pos), __impl(impl) { }
    __pfor_vector_iterator()
      : __pos(), __impl() { }
    
    value_type operator*() const { return __impl->operator[](__pos); }
    
    self_type& operator++() { ++ __pos; return *this; }
    self_type& operator--() { -- __pos; return *this; }
    
    self_type& operator+=(difference_type __n) { __pos += __n; return *this; }
    self_type& operator-=(difference_type __n) { __pos -= __n; return *this; }
    
    self_type operator++(int) { self_type __tmp = *this; ++ *this; return __tmp; }
    self_type operator--(int) { self_type __tmp = *this; -- *this; return __tmp; }
    self_type operator+(difference_type __n) const { self_type __tmp = *this; return __tmp += __n; }
    self_type operator-(difference_type __n) const { self_type __tmp = *this; return __tmp -= __n; }
    
  private:
    size_type        __pos;
    const impl_type* __impl;
  };
  
  template <typename Tp, typename Impl>
  inline bool
  operator==(const __pfor_vector_iterator<Tp,Impl>& x,
	     const __pfor_vector_iterator<Tp,Impl>& y)
  {
    return x.__impl == y.__impl && x.__pos == y.__pos;
  }
  
  template <typename Tp, typename Impl>
  inline bool
  operator!=(const __pfor_vector_iterator<Tp,Impl>& x,
	     const __pfor_vector_iterator<Tp,Impl>& y)
  {
    return x.__impl != y.__impl || x.__pos != y.__pos;
  }



  template <typename Tp, typename Iterator>
  class __pfor_pseudo_vector
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef int64_t   off_type;
    typedef uint32_t  pos_type;
    typedef Tp        value_type;
    
  public:
    __pfor_pseudo_vector() {}
    __pfor_pseudo_vector(Iterator first, Iterator last)
      : __first(first), __last(last) {}
    
    size_type size() const { return __last - __first; }
    bool empty() const { return __first == __last; }
    
    Iterator begin() const { return __first; }
    Iterator end() const { return __last; }
    
    const value_type& operator[](size_type pos) const { return *(__first + pos); }
    
  private:
    Iterator __first;
    Iterator __last;
  };
  

  template <typename Tp>
  struct __pfor_vector_base
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef int64_t   off_type;
    typedef uint32_t  pos_type;
    
    typedef Tp        value_type;

    typedef boost::filesystem::path path_type;

    static const size_type entry_size = 128;         // 128 values... must be 128!
    static const size_type block_size = 1024 * 1024; // 1024K values... must be power of two
    
    static const size_type entry_bits = utils::bithack::static_floor_log2<entry_size>::result;
    static const size_type block_bits = utils::bithack::static_floor_log2<block_size>::result;
    
    static const size_type entry_mask = (1 << entry_bits) - 1;
    static const size_type block_mask = (1 << block_bits) - 1;

    static const size_type entry_per_block = block_size >> entry_bits;
    
    struct header_type
    {
      value_type base;
      pos_type   bits;
      pos_type   pos_code;
      pos_type   pos_exception;
      
      header_type() {}
      
      header_type(const value_type& __base,
		  const pos_type& __bits,
		  const pos_type& __pos_code,
		  const pos_type& __pos_exception)
	: base(__base), bits(__bits), pos_code(__pos_code), pos_exception(__pos_exception) {}
    };
    
    typedef uint32_t   entry_type;
  };
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class pfor_vector_mapped : public __pfor_vector_base<Tp>
  {
  private:
    typedef __pfor_vector_base<Tp> base_type;
  public:
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::off_type        off_type;
    typedef typename base_type::pos_type        pos_type;
    typedef typename base_type::path_type       path_type;
    typedef typename base_type::value_type      value_type;
    
  private:
    
    using base_type::entry_size;
    using base_type::block_size;
    using base_type::entry_bits;
    using base_type::block_bits;
    using base_type::entry_mask;
    using base_type::block_mask;
    using base_type::entry_per_block;

    typedef typename base_type::header_type header_type;
    typedef typename base_type::entry_type  entry_type;
    
    typedef Alloc value_alloc_type;
    typedef typename Alloc::template rebind<entry_type>::other  entry_alloc_type;
    typedef typename Alloc::template rebind<header_type>::other header_alloc_type;
    
    typedef utils::map_file<header_type, header_alloc_type> header_set_type;
    typedef utils::map_file<entry_type, entry_alloc_type>   entry_set_type;
    typedef utils::map_file<value_type, value_alloc_type>   code_set_type;
    typedef utils::map_file<value_type, value_alloc_type>   exception_set_type;
    
    typedef codec::pfor_codec<value_type> codec_type;
    typedef __pfor_value_proxy<value_type> proxy_type;
    
    typedef __pfor_pseudo_vector<entry_type, typename entry_set_type::const_iterator>     pseudo_entry_set_type;
    typedef __pfor_pseudo_vector<value_type, typename code_set_type::const_iterator>      pseudo_code_set_type;
    typedef __pfor_pseudo_vector<value_type, typename exception_set_type::const_iterator> pseudo_exception_set_type;

  private:
    typedef pfor_vector_mapped<Tp, Alloc> self_type;
    
  public:
    typedef __pfor_vector_iterator<Tp, self_type>       iterator;
    typedef __pfor_vector_iterator<Tp, self_type> const_iterator;
    
  public:
    pfor_vector_mapped() {}
    pfor_vector_mapped(const path_type& path) { open(path); }
    
  public:
    size_type size() const { return __size; }
    bool empty() const { return __size == 0; }
    path_type path() const { return __headers.path().parent_path(); }
    
    const_iterator begin() const { return const_iterator(0, this); }
    const_iterator end() const { return const_iterator(size(), this); }
    
    value_type operator[](size_type pos) const
    {
      const size_type header_pos = pos >> block_bits;
      const size_type block_pos = pos & block_mask;
      
      const pos_type pos_code = __headers[header_pos].pos_code;
      const pos_type pos_exception = __headers[header_pos].pos_exception;
      
      return codec_type::decompress(block_pos,
				    __headers[header_pos].base,
				    __headers[header_pos].bits,
				    pseudo_entry_set_type(__entries.begin() + header_pos * entry_per_block, __entries.end()),
				    pseudo_code_set_type(__codes.begin() + pos_code, __codes.end()),
				    pseudo_exception_set_type(__exceptions.begin() + pos_exception, __exceptions.end()));
    }
    
    void clear()
    {
      __headers.clear();
      __entries.clear();
      __codes.clear();
      __exceptions.clear();
      __size = 0;
    }
    void close() { clear(); }
    
    void open(const path_type& path)
    {
      typedef utils::repository repository_type;
      
      clear();
      
      repository_type rep(path, repository_type::read);
      
      __headers.open(rep.path("header"));
      __entries.open(rep.path("entry"));
      __codes.open(rep.path("code"));
      __exceptions.open(rep.path("exception"));
      
      repository_type::const_iterator iter = rep.find("size");
      if (iter == rep.end())
	throw std::runtime_error("no size?");
      __size = atoll(iter->second.c_str());
    }

    
  public:
    header_set_type    __headers;
    entry_set_type     __entries;
    code_set_type      __codes;
    exception_set_type __exceptions;
    size_type          __size;
  };

  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class pfor_vector : public __pfor_vector_base<Tp>
  {
  private:
    typedef __pfor_vector_base<Tp> base_type;
    
  public:
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::off_type        off_type;
    typedef typename base_type::pos_type        pos_type;
    typedef typename base_type::path_type       path_type;
    typedef typename base_type::value_type      value_type;
    
  private:
    using base_type::entry_size;
    using base_type::block_size;
    using base_type::entry_bits;
    using base_type::block_bits;
    using base_type::entry_mask;
    using base_type::block_mask;
    using base_type::entry_per_block;

    typedef typename base_type::header_type header_type;
    typedef typename base_type::entry_type  entry_type;
    
    typedef Alloc value_alloc_type;
    typedef typename Alloc::template rebind<entry_type>::other  entry_alloc_type;
    typedef typename Alloc::template rebind<header_type>::other header_alloc_type;
    
    typedef std::vector<value_type, value_alloc_type> value_set_type;

    typedef std::vector<header_type, header_alloc_type> header_set_type;
    typedef std::vector<entry_type, entry_alloc_type>   entry_set_type;
    typedef std::vector<value_type, value_alloc_type>   code_set_type;
    typedef std::vector<value_type, value_alloc_type>   exception_set_type;
    
    typedef codec::pfor_codec<value_type> codec_type;
    typedef __pfor_value_proxy<value_type> proxy_type;
    
    typedef __pfor_pseudo_vector<entry_type, typename entry_set_type::const_iterator>     pseudo_entry_set_type;
    typedef __pfor_pseudo_vector<value_type, typename code_set_type::const_iterator>      pseudo_code_set_type;
    typedef __pfor_pseudo_vector<value_type, typename exception_set_type::const_iterator> pseudo_exception_set_type;

    
  private:
    typedef pfor_vector<Tp, Alloc> self_type;
    
  public:
    typedef __pfor_vector_iterator<Tp, self_type>       iterator;
    typedef __pfor_vector_iterator<Tp, self_type> const_iterator;
    
  public:
    pfor_vector() { clear(); }
    
  public:
    size_type size() const { return (! __values.empty() ? __values.size() : __size); }
    bool empty() const { return __values.empty() && __size == 0; }
    
    const_iterator begin() const { return const_iterator(0, this); }
    const_iterator end() const { return const_iterator(size(), this); }
    
    value_type operator[](size_type pos) const
    {
      if (! __values.empty()) {
	return __values[pos];
      } else {
	
	const size_type header_pos = pos >> block_bits;
	const size_type block_pos = pos & block_mask;
	
	const pos_type pos_code = __headers[header_pos].pos_code;
	const pos_type pos_exception = __headers[header_pos].pos_exception;
	
	return codec_type::decompress(block_pos,
				      __headers[header_pos].base,
				      __headers[header_pos].bits,
				      pseudo_entry_set_type(__entries.begin() + header_pos * entry_per_block, __entries.end()),
				      pseudo_code_set_type(__codes.begin() + pos_code, __codes.end()),
				      pseudo_exception_set_type(__exceptions.begin() + pos_exception, __exceptions.end()));
      }
    }


    void push_back(const value_type& x)
    {
      if (__size > 0)
	throw std::runtime_error("you cannot add data to the compressed block");
      
      __values.push_back(x);
    }
    
    template <typename Iterator>
    void insert(Iterator first, Iterator last)
    {
      if (__size > 0)
	throw std::runtime_error("you cannot add data to the compressed block");
      
      __values.insert(__values.end(), first, last);
    }
    
    void clear()
    {
      __values.clear();
      
      __headers.clear();
      __entries.clear();
      __codes.clear();
      __exceptions.clear();
      __size = 0;
    }
    void close() { clear(); }
    
    
    void build()
    {
      __headers.clear();
      __entries.clear();
      __codes.clear();
      __exceptions.clear();
      
      __size = __values.size();
      
      if (__values.empty())
	return;

      // we will prepare local entries, codes, exceptions...
      entry_set_type     entries;
      code_set_type      codes;
      exception_set_type exceptions;
      
      const size_type num_block = __values.size() >> block_bits;
      for (int i = 0; i < num_block; ++ i) {
	
	std::pair<value_type, size_type> result = codec_type::estimate(&(*(__values.begin() + i * block_size)),
								       block_size);
	
	codec_type::compress(&(*(__values.begin() + i * block_size)),
			     block_size,
			     result.first,
			     result.second,
			     entries,
			     codes,
			     exceptions);
	

	__headers.push_back(header_type(result.first, result.second, __codes.size(), __exceptions.size()));
	__entries.insert(__entries.end(), entries.begin(), entries.end());
	__codes.insert(__codes.end(), codes.begin(), codes.end());
	__exceptions.insert(__exceptions.end(), exceptions.begin(), exceptions.end());
	
	entries.clear();
	codes.clear();
	exceptions.clear();
      }
      
      if (__values.size() & block_mask) {
	// compress the last block...
	std::pair<value_type, size_type> result = codec_type::estimate(&(*(__values.begin() + num_block * block_size)),
								       __values.size() & block_mask);
	
	codec_type::compress(&(*(__values.begin() + num_block * block_size)),
			     __values.size() & block_mask,
			     result.first,
			     result.second,
			     entries,
			     codes,
			     exceptions);
	
	__headers.push_back(header_type(result.first, result.second, __codes.size(), __exceptions.size()));
	__entries.insert(__entries.end(), entries.begin(), entries.end());
	__codes.insert(__codes.end(), codes.begin(), codes.end());
	__exceptions.insert(__exceptions.end(), exceptions.begin(), exceptions.end());
	
	entries.clear();
	codes.clear();
	exceptions.clear();
      }
	
      __values.clear();
    }
    
    void write(const path_type& path)
    {
      typedef utils::repository repository_type;
      
      if (! __values.empty())
	build();
      
      repository_type rep(path, repository_type::write);
      
      dump_file(rep.path("header"), __headers);
      dump_file(rep.path("entry"), __entries);
      dump_file(rep.path("code"), __codes);
      dump_file(rep.path("exception"), __exceptions);
      
      std::ostringstream stream_size;
      std::ostringstream stream_integral_size;
      std::ostringstream stream_compressed_size;
      stream_size << __size;
      stream_integral_size << sizeof(Tp);
      stream_compressed_size << (__headers.size() * sizeof(header_type)
				 + __entries.size() * sizeof(entry_type)
				 + __codes.size() * sizeof(value_type)
				 + __exceptions.size() * sizeof(value_type));
      rep["size"] = stream_size.str();
      rep["integral-size"] = stream_integral_size.str();
      rep["compressed-size"] = stream_compressed_size.str();
      rep["type"] = "pfor";
    }
    
  private:
    template <typename _Path, typename _Data>
    inline
    void dump_file(const _Path& file, const _Data& data)
    {
      std::auto_ptr<boost::iostreams::filtering_ostream> os(new boost::iostreams::filtering_ostream());
      os->push(boost::iostreams::file_sink(file.native_file_string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
      
      const int64_t file_size = sizeof(typename _Data::value_type) * data.size();
      for (int64_t offset = 0; offset < file_size; offset += 1024 * 1024)
	os->write(((char*) &(*data.begin())) + offset, std::min(int64_t(1024 * 1024), file_size - offset));
    }

    
  public:
    value_set_type     __values;

    header_set_type    __headers;
    entry_set_type     __entries;
    code_set_type      __codes;
    exception_set_type __exceptions;
    size_type          __size;
  };
};

#endif
