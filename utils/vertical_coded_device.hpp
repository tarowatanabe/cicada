// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__VERTICAL_CODED_DEVICE__HPP__
#define __UTILS__VERTICAL_CODED_DEVICE__HPP__ 1

#include <vector>

#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/filesystem/path.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/file.hpp>

#include <utils/repository.hpp>

#include <boost/shared_ptr.hpp>

namespace utils
{
  
  // vertical coded sink
  // we assume fixed size insertion...
  // when closed, performs building of compressed stream..  
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class vertical_coded_sink
  {
  public:
    typedef size_t    size_type;
    typedef char      char_type;
    typedef uint8_t   byte_type;
    typedef uint64_t  off_type;
    
    typedef Tp        value_type;
    
    typedef boost::filesystem::path path_type;
    
    struct category : public boost::iostreams::sink_tag,
		      public boost::iostreams::closable_tag {};
    
    vertical_coded_sink(const path_type& file, const size_type buffer_size=4096) : pimpl(new impl()) { open(file, buffer_size); }
    vertical_coded_sink() : pimpl(new impl()) {}
    
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
    void open(const path_type& file, const size_type buffer_size=4096) { pimpl->open(file, buffer_size); }
    void close() { pimpl->close(); }
    
  private:
    struct impl
    {
      typedef typename Alloc::template rebind<byte_type>::other  byte_alloc_type;
      typedef typename Alloc::template rebind<value_type>::other value_alloc_type;
      
      typedef std::vector<byte_type,  byte_alloc_type>  buffer_type;
      typedef std::vector<value_type, value_alloc_type> inverted_type;
      
      typedef boost::iostreams::filtering_ostream ostream_type;
      
      impl() {}
      ~impl() { close(); }
      
      std::streamsize write(const char_type* s, std::streamsize n);
      void open(const path_type& file, const size_type buffer_size=4096);
      bool is_open() const;
      void close();
      
      buffer_type buffer;
      size_type   offset;
      off_type   size;
      
      boost::shared_ptr<ostream_type> os_data;
      boost::shared_ptr<ostream_type> os_off;
      inverted_type inverted;
      path_type path;
    };
    
    boost::shared_ptr<impl> pimpl;
  };
  
  template <typename Tp,typename Alloc>
  inline
  bool vertical_coded_sink<Tp,Alloc>::impl::is_open() const { return os_data; }
  
  template <typename Tp,typename Alloc>
  inline
  std::streamsize vertical_coded_sink<Tp,Alloc>::impl::write(const char_type* s, std::streamsize n)
  {
    if (! is_open()) return -1;
    
    const std::streamsize copy_size = std::min(std::streamsize(buffer.size() - offset), n);
    
    std::copy(s, s + copy_size, buffer.begin() + offset);
    offset += copy_size;

    if (offset == buffer.size()) {
      const value_type* first = reinterpret_cast<const value_type*>(&(*buffer.begin()));
      const value_type* last = reinterpret_cast<const value_type*>(&(*buffer.end()));
      
      for (const value_type* iter = first; iter != last; ++ iter, ++ size) {
	const byte_type  value  = (*iter) & 0xff;
	const value_type invert = (*iter) >> 8;
	
	os_data->write((char*) &value, sizeof(byte_type));
	
	if (inverted.empty() || value_type(inverted.size() - 1) != invert)
	  inverted.resize(invert + 1, size);
      }
      
      offset = 0;
    }
    
    return copy_size;
  }
  
  template <typename Tp,typename Alloc>
  inline
  void vertical_coded_sink<Tp,Alloc>::impl::close()
  {
    typedef utils::repository repository_type;

    if (is_open()) {
      if (offset > 0) {
	if (offset % sizeof(value_type) != 0)
	  throw std::runtime_error("invalid write to this device...");
	
	buffer.resize(offset);
	
	const value_type* first = reinterpret_cast<const value_type*>(&(*buffer.begin()));
	const value_type* last = reinterpret_cast<const value_type*>(&(*buffer.end()));
	
	for (const value_type* iter = first; iter != last; ++ iter, ++ size) {
	  const byte_type  value  = (*iter) & 0xff;
	  const value_type invert = (*iter) >> 8;
	  
	  os_data->write((char*) &value, sizeof(byte_type));
	  
	  if (inverted.empty() || value_type(inverted.size() - 1) != invert)
	    inverted.resize(invert + 1, size);
	}
      }
      
      repository_type rep(path, repository_type::read);
      os_off.reset(new boost::iostreams::filtering_ostream());
      os_off->push(boost::iostreams::file_sink(rep.path("offsets").file_string()));
      os_off->exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
      
      off_type off(0);
      os_off->write((char*) &off, sizeof(off_type));
      os_off->write((char*) &size, sizeof(off_type));

      const off_type value_size = size;

      bool finished = (inverted.back() <= 0xff);
      int code_size = 1;
      
      inverted_type inverted_new;
      while (! finished) {
	for (size_type i = 0; i < inverted.size(); ++ i, ++ size) {
	  const byte_type  value  = inverted[i] & 0xff;
	  const value_type invert = inverted[i] >> 8;
	  
	  os_data->write((char*) &value, sizeof(byte_type));
	  
	  if (inverted_new.empty() || value_type(inverted_new.size() - 1) != invert)
	    inverted_new.resize(invert + 1, i);
	}
	
	inverted_new.swap(inverted);
	inverted_new.clear();
	
	os_off->write((char*) &size, sizeof(off_type));
	++ code_size;
	
	if (inverted.back() <= 0xff)
	  finished = true;
      }
      
      for (size_type i = 0; i < inverted.size(); ++ i, ++ size) {
	const byte_type value = inverted[i] & 0xff;
	os_data->write((char*) &value, sizeof(byte_type));
      }
      os_off->write((char*) &size, sizeof(off_type));
      ++ code_size;
      
      
      // repository...
      {	
	std::ostringstream stream_integral_size;
	std::ostringstream stream_size;
	std::ostringstream stream_coded_size;
	
	stream_integral_size << sizeof(value_type);
	stream_size << value_size;
	stream_coded_size << code_size;
	
	rep["size"] = stream_size.str();
	rep["integral-size"] = stream_integral_size.str();
	rep["code-size"] = stream_coded_size.str();
      }
    }
    
    buffer.clear();
    offset = 0;
    size = 0;
    os_data.reset();
    os_off.reset();
    inverted.clear();
    path = path_type();
  }
  
  template <typename Tp,typename Alloc>
  inline
  void vertical_coded_sink<Tp,Alloc>::impl::open(const path_type& file, const size_type buffer_size)
  {
    typedef utils::repository repository_type;
    
    close();
      
    repository_type rep(file, repository_type::write);
    rep["type"] = "vertical-coded";
    
    os_data.reset(new boost::iostreams::filtering_ostream());
    os_data->push(boost::iostreams::file_sink(rep.path("data").file_string()), buffer_size);
    os_data->exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
    
    const size_type buffer_size_new = ((buffer_size + sizeof(value_type) - 1) / sizeof(value_type)) * sizeof(value_type);
    
    buffer.reserve(buffer_size_new);
    buffer.resize(buffer_size_new);
    offset = 0;
    size = 0;
    inverted.clear();
    path = file;
  }
};

#endif
