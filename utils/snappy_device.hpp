// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SNAPPY_DEVICE__HPP__
#define __UTILS__SNAPPY_DEVICE__HPP__ 1

#include <vector>

#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/filesystem/path.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/file.hpp>

#include <utils/vertical_coded_device.hpp>
#include <utils/repository.hpp>

#include <boost/shared_ptr.hpp>

#include <snappy.hpp>

namespace utils
{
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class snappy_sink
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
    
    snappy_sink(const path_tyep& file) : pimpl(new impl()) { open(file); }
    snappy_sink() : pimpl(new impl()) {}

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
      typedef typename Alloc::template rebind<byte_type>::other  byte_alloc_type;
      typedef std::vector<byte_type,  byte_alloc_type>  buffer_type;
      typedef boost::iostreams::filtering_ostream ostream_type;
      
      impl() {}
      ~impl() { close(); }

      std::streamsize write(const char_type* s, std::streamsize n);
      void open(const path_type& file);
      bool is_open() const;
      void close();
      
      buffer_type buffer;
      buffer_type buffer_compressed;
      
      off_type offset;
      off_type offset_buffer;
      off_type offset_compressed;
      
      boost::shared_ptr<ostream_type> os_data;
      boost::shared_ptr<ostream_type> os_offset;
      
      path_type path;
    };
    
    boost::shared_ptr<impl> pimpl;
  };

  template <typename Tp,typename Alloc>
  inline
  bool snappy_sink<Tp,Alloc>::impl::is_open() const { return os_data; }

  template <typename Tp, typename Alloc>
  inline
  std::streamsize snappy_sink<Tp,Alloc>::impl::write(const char_type* s, std::streamsize n)
  {
    if (! is_open()) return -1;
    
    const std::streamsize copy_size = std::min(std::streamsize(buffer.size() - offset_buffer), n);
    
    std::copy(s, s + copy_size, buffer.begin() + offset_buffer);
    offset_buffer += copy_size;
    offset += copy_size;
    
    if (offset_buffer == buffer.size()) {
      // compress and dump in os_data
      
      size_t compressed_length = 0;
      snappy::RawCompress(&(*buffer.begin()), offset_buffer, &(*buffer_compressed.begin()), compressed_length);
      
      offset_compressed += compressed_length;
      offset_buffer = 0;
      
      os_data->write((char*) &(*buffer_compressed.begin()), compressed_length);
      os_offset->write((char*) &offset_compressed, sizeof(offset_compressed));
    }
    
    return copy_size;
  }
  
  template <typename Tp,typename Alloc>
  inline
  void snappy_sink<Tp,Alloc>::impl::close()
  {
    typedef utils::repository repository_type;
    
    if (is_open()) {
      if (offset_buffer) {
	size_t compressed_length = 0;
	snappy::RawCompress(&(*buffer.begin()), offset_buffer, &(*buffer_compressed.begin()), compressed_length);
	
	offset_compressed += compressed_length;
	offset_buffer = 0;
	
	os_data->write((char*) &(*buffer_compressed.begin()), compressed_length);
	os_offset->write((char*) &offset_compressed, sizeof(offset_compressed));
      }
      
      repository_type rep(path, repository_type::read);
      
      std::ostringstream stream_size;
      std::ostringstream stream_buffer;
      std::ostringstream stream_integral_size;
      stream_size << offset;
      stream_buffer << buffer.size();
      stream_integral_size << sizeof(Tp);
      
      rep["size"] = stream_size.str();
      rep["block"] = stream_buffer.str();
      rep["integral-size"] = stream_integral_size.str();
    }
    
    buffer.clear();
    buffer_compressed.clear();
    
    offset = 0;
    offset_buffer = 0;
    offset_compressed = 0;
    
    os_data.reset();
    os_offset.reset();

    path = path_type();
  }

  template <typename Tp,typename Alloc>
  inline
  void snappy_sink<Tp,Alloc>::impl::open(const path_type& file)
  {
    typedef utils::repository repository_type;
    typedef typename Alloc::template rebind<off_type>::other  off_alloc_type;
    
    close();
    
    repository_type rep(file, repository_type::write);
    rep["type"] = "snappy";
    
    const size_type max_compressed_length = snappy::MaxCompressedLength(4096);
    
    buffer.reserve(4096);
    buffer.resize(4096);
    
    buffer_compressed.reserve(max_compressed_length);
    buffer_compressed.resize(max_compressed_length);
    
    offset = 0;
    offset_buffer = 0;
    offset_compressed = 0;

    os_data.reset(new ostream_type());
    os_data->push(boost::iostreams::file_sink(rep.path("data").string()), 1024 * 1024);
    os_data->exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
    
    os_offset.reset(new ostream_type());
    os_offset->push(utils::vertical_coded_sink<off_type, off_alloc_type>(rep.path("offsets")));
    os_offset->exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
    
    os_offset->write((char*) &offset_compressed, sizeof(offset_compressed));
        
    path = file;
  }
};

#endif
