// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CODEC__SNAPPY_HPP__
#define __CODEC__SNAPPY_HPP__ 1

#include <cstddef>
#include <memory>

#include <boost/iostreams/filter/symmetric.hpp>

namespace codec
{
  namespace detail
  {
    struct snappy_param
    {
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      typedef char byte_type;
      typedef char char_type;

      static const size_type chunk_size;
      static const size_type bound_size;
    };
  };

  struct snappy_compressor_impl : detail::snappy_param
  {
    //
    // we will initialize buffer, but do not allow copy/assign
    //
    snappy_compressor_impl();
    snappy_compressor_impl(const snappy_compressor_impl&);
    ~snappy_compressor_impl()
    {
      delete [] buffer;
      delete [] buffer_compressed;
    }
    
    snappy_compressor_impl& operator=(const snappy_compressor_impl&) { return *this; }
    
    bool filter(const char*& src_begin,
		const char* src_end,
                char*& dest_begin,
		char* dest_end,
		bool flush);
    
    void close()
    {
      pos = 0;
      size_compressed = 0;
      pos_compressed = 0;
    }
    
    byte_type* buffer;
    size_type  pos;
    
    byte_type* buffer_compressed;
    size_type  size_compressed;
    size_type  pos_compressed;
  };
  
  struct snappy_decompressor_impl : detail::snappy_param
  {
    snappy_decompressor_impl();
    snappy_decompressor_impl(const snappy_decompressor_impl&);
    ~snappy_decompressor_impl()
    {
      delete [] buffer_compressed;
      delete [] buffer;
    }
    
    snappy_decompressor_impl& operator=(const snappy_decompressor_impl&) { return *this; }
    
    bool filter(const char*& begin_in,
		const char* end_in,
		char*& begin_out,
		char* end_out,
		bool flush);
    void close()
    {
      size_compressed = 0;
      pos_compressed = 0;
      size = 0;
      pos = 0;
    }
    
    byte_type* buffer_compressed;
    size_type  size_compressed;
    size_type  pos_compressed;
    
    byte_type* buffer;
    size_type  size;
    size_type  pos;
  };

  template<typename Alloc = std::allocator<char> >
  struct basic_snappy_compressor
    : boost::iostreams::symmetric_filter<snappy_compressor_impl, Alloc>
  {
  private:
    typedef snappy_compressor_impl                               impl_type;
    typedef boost::iostreams::symmetric_filter<impl_type, Alloc> base_type;
    
  public:
    typedef typename base_type::char_type char_type;
    typedef typename base_type::category  category;

  public:
    basic_snappy_compressor(int buffer_size = boost::iostreams::default_device_buffer_size)
      : base_type(buffer_size) {}
  };
  
  template<typename Alloc = std::allocator<char> >
  struct basic_snappy_decompressor
    : boost::iostreams::symmetric_filter<snappy_decompressor_impl, Alloc>
  {
  private:
    typedef snappy_decompressor_impl                             impl_type;
    typedef boost::iostreams::symmetric_filter<impl_type, Alloc> base_type;
    
  public:
    typedef typename base_type::char_type char_type;
    typedef typename base_type::category  category;

  public:
    basic_snappy_decompressor(int buffer_size = boost::iostreams::default_device_buffer_size)
      : base_type(buffer_size) {}
  };

  typedef basic_snappy_compressor<>   snappy_compressor;
  typedef basic_snappy_decompressor<> snappy_decompressor;
};

#endif
