// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CODEC__QUICKLZ_HPP__
#define __CODEC__QUICKLZ_HPP__ 1

#include <cstddef>
#include <memory>

#include <boost/iostreams/filter/symmetric.hpp>

namespace codec
{
  namespace detail
  {
    struct quicklz_param_impl;

    struct quicklz_param
    {
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      typedef char byte_type;
      typedef char char_type;
      
      static const size_type chunk_size;
      static const size_type bound_size;
      
      quicklz_param();
      quicklz_param(const quicklz_param&x);
      ~quicklz_param();
      quicklz_param& operator=(const quicklz_param&) { return *this; }
      
      quicklz_param_impl* pimpl;
    };
  };

  struct quicklz_compressor_impl : detail::quicklz_param
  {
    //
    // we will initialize buffer, but do not allow copy/assign
    //
    quicklz_compressor_impl();
    quicklz_compressor_impl(const quicklz_compressor_impl&);
    ~quicklz_compressor_impl()
    {
      delete [] buffer;
      delete [] buffer_compressed;
    }
    
    quicklz_compressor_impl& operator=(const quicklz_compressor_impl&) { return *this; }
    
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
  
  struct quicklz_decompressor_impl : detail::quicklz_param
  {
    quicklz_decompressor_impl();
    quicklz_decompressor_impl(const quicklz_decompressor_impl&);
    ~quicklz_decompressor_impl()
    {
      delete [] buffer_compressed;
      delete [] buffer;
    }
    
    quicklz_decompressor_impl& operator=(const quicklz_decompressor_impl&) { return *this; }
    
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
  struct basic_quicklz_compressor
    : boost::iostreams::symmetric_filter<quicklz_compressor_impl, Alloc>
  {
  private:
    typedef quicklz_compressor_impl                              impl_type;
    typedef boost::iostreams::symmetric_filter<impl_type, Alloc> base_type;
    
  public:
    typedef typename base_type::char_type char_type;
    typedef typename base_type::category  category;

  public:
    basic_quicklz_compressor(int buffer_size = boost::iostreams::default_device_buffer_size)
      : base_type(buffer_size) {}
  };
  
  template<typename Alloc = std::allocator<char> >
  struct basic_quicklz_decompressor
    : boost::iostreams::symmetric_filter<quicklz_decompressor_impl, Alloc>
  {
  private:
    typedef quicklz_decompressor_impl                            impl_type;
    typedef boost::iostreams::symmetric_filter<impl_type, Alloc> base_type;
    
  public:
    typedef typename base_type::char_type char_type;
    typedef typename base_type::category  category;

  public:
    basic_quicklz_decompressor(int buffer_size = boost::iostreams::default_device_buffer_size)
      : base_type(buffer_size) {}
  };

  typedef basic_quicklz_compressor<>   quicklz_compressor;
  typedef basic_quicklz_decompressor<> quicklz_decompressor;
};

#endif
