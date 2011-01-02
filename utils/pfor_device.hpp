// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// a device implemenation for pfor-vector

#ifndef __UTILS__PFOR_DEVICE__HPP__
#define __UTILS__PFOR_DEVICE__HPP__ 1

#include <stdint.h>

#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>

#include <boost/filesystem.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/file.hpp>

#include <utils/repository.hpp>
#include <utils/bithack.hpp>
#include <utils/bitpack.hpp>
#include <utils/pfor_vector.hpp>

#include <codec/pfor_codec.hpp>

namespace utils
{

  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class pfor_sink : public __pfor_vector_base<Tp>
  {
  public:
    typedef size_t    size_type;
    typedef char      char_type;
    typedef uint8_t   byte_type;
    
    struct category : public boost::iostreams::sink_tag,
		      public boost::iostreams::closable_tag {};
    
  private:
    typedef __pfor_vector_base<Tp> base_type;
    
  public:
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
    
  public:
    pfor_sink(const path_type& file) : pimpl(new impl()) { open(file); }
    pfor_sink() : pimpl(new impl()) {}
    
  public:
    std::streamsize write(const char_type* s, std::streamsize n)
    {
      std::streamsize __n = 0;
      std::streamsize available = n;
      while (available) {
	const std::streamsize written = pimpl->write(s, available);
	if (written == -1)
	  return (__n > 0 ? __n : -1);
	available -= written;
	s += written;
	__n += written;
      }
      return __n;
    }
    
    bool is_open() const { return pimpl->is_open(); }
    void open(const path_type& file) { pimpl->open(file); }
    void close() { pimpl->close(); }
    
  private:
    struct impl
    {
      typedef boost::iostreams::filtering_ostream ostream_type;
      
      typedef typename Alloc::template rebind<byte_type>::other  byte_alloc_type;
      typedef std::vector<byte_type,  byte_alloc_type>  buffer_type;
      
      impl() {}
      ~impl() { close(); }
      
      std::streamsize write(const char_type* s, std::streamsize n);
      void open(const path_type& file);
      bool is_open() const;
      void close();
      
      buffer_type buffer;
      size_type offset;
      size_type size;
      
      boost::shared_ptr<ostream_type> os_header;
      boost::shared_ptr<ostream_type> os_entry;
      boost::shared_ptr<ostream_type> os_code;
      boost::shared_ptr<ostream_type> os_exception;
      
      pos_type code_size;
      pos_type exception_size;
      
      entry_set_type     entries;
      code_set_type      codes;
      exception_set_type exceptions;

      path_type path;
    };
    
    boost::shared_ptr<impl> pimpl;
  };
  
  
  template <typename Tp,typename Alloc>
  inline
  bool pfor_sink<Tp,Alloc>::impl::is_open() const
  {
    return os_header;
  } 
  
  template <typename Tp,typename Alloc>
  inline
  std::streamsize pfor_sink<Tp,Alloc>::impl::write(const char_type* s, std::streamsize n)
  {
    if (! is_open()) return -1;
    
    const std::streamsize copy_size = std::min(std::streamsize(buffer.size() - offset), n);
    
    std::copy(s, s + copy_size, buffer.begin() + offset);
    offset += copy_size;
    
    if (offset == buffer.size()) {
      entries.clear();
      codes.clear();
      exceptions.clear();

      std::pair<value_type, size_type> result = codec_type::estimate(reinterpret_cast<const value_type*>(&(*buffer.begin())),
								     block_size);
      codec_type::compress(reinterpret_cast<const value_type*>(&(*buffer.begin())),
			   block_size,
			   result.first,
			   result.second,
			   entries,
			   codes, exceptions);
      
      header_type header(result.first, result.second, code_size, exception_size);
      
      os_header->write((char*) &header, sizeof(header_type));
      os_entry->write((char*) &(*entries.begin()), sizeof(entry_type) * entries.size());
      os_code->write((char*) &(*codes.begin()), sizeof(value_type) * codes.size());
      os_exception->write((char*) &(*exceptions.begin()), sizeof(value_type) * exceptions.size());
      
      code_size += codes.size();
      exception_size += exceptions.size();
      size += block_size;
      
      entries.clear();
      codes.clear();
      exceptions.clear();
      
      offset = 0;
    }
    
    return copy_size;
  }

  template <typename Tp,typename Alloc>
  inline
  void pfor_sink<Tp,Alloc>::impl::close()
  {
    if (is_open()) {
      
      if (offset > 0) {
	if (offset % sizeof(value_type) != 0)
	  throw std::runtime_error("invalid write to this device...");
	
	entries.clear();
	codes.clear();
	exceptions.clear();

	std::pair<value_type, size_type> result = codec_type::estimate(reinterpret_cast<const value_type*>(&(*buffer.begin())),
								       offset / sizeof(value_type));
	codec_type::compress(reinterpret_cast<const value_type*>(&(*buffer.begin())),
			     offset / sizeof(value_type),
			     result.first,
			     result.second,
			     entries,
			     codes, exceptions);
      
	header_type header(result.first, result.second, code_size, exception_size);
	
	os_header->write((char*) &header, sizeof(header_type));
	os_entry->write((char*) &(*entries.begin()), sizeof(entry_type) * entries.size());
	os_code->write((char*) &(*codes.begin()), sizeof(value_type) * codes.size());
	os_exception->write((char*) &(*exceptions.begin()), sizeof(value_type) * exceptions.size());
	
	code_size += codes.size();
	exception_size += exceptions.size();
	size += offset / sizeof(value_type);
	
	entries.clear();
	codes.clear();
	exceptions.clear();
	
	offset = 0;
      }
      
      {
	utils::repository rep(path, utils::repository::read);
	std::ostringstream stream_size;
	std::ostringstream stream_integral_size;
	stream_size << size;
	stream_integral_size << sizeof(Tp);
	rep["size"] = stream_size.str();
	rep["integral-size"] = stream_integral_size.str();
      }
    }
    
    buffer.clear();
    offset = 0;
    size = 0;
    
    os_header.reset();
    os_entry.reset();
    os_code.reset();
    os_exception.reset();
    
    code_size = 0;
    exception_size = 0;
    
    entries.clear();
    codes.clear();
    exceptions.clear();
  }
  
  template <typename Tp,typename Alloc>
  inline
  void pfor_sink<Tp,Alloc>::impl::open(const path_type& file)
  {
    close();
    
    buffer.reserve(block_size * sizeof(value_type));
    buffer.resize(block_size * sizeof(value_type));
    
    utils::repository rep(file, utils::repository::write);
    rep["type"] = "pfor";
    
    os_header.reset(new ostream_type());
    os_entry.reset(new ostream_type());
    os_code.reset(new ostream_type());
    os_exception.reset(new ostream_type());
    
    os_header->push(boost::iostreams::file_sink(rep.path("header").file_string()));
    os_entry->push(boost::iostreams::file_sink(rep.path("entry").file_string()));
    os_code->push(boost::iostreams::file_sink(rep.path("code").file_string()));
    os_exception->push(boost::iostreams::file_sink(rep.path("exception").file_string()));

    path = file;    
  }
};

#endif
