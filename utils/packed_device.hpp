// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__PACKED_DEVICE__HPP__
#define __UTILS__PACKED_DEVICE__HPP__ 1

#include <vector>

#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/filesystem/path.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/file.hpp>

#include <utils/succinct_vector.hpp>
#include <utils/packed_vector.hpp>
#include <utils/repository.hpp>

#include <boost/shared_ptr.hpp>

namespace utils
{
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class packed_sink
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
    
    packed_sink(const path_type& file) : pimpl(new impl()) { open(file); }
    packed_sink() : pimpl(new impl()) {}
    
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
    struct impl : public __packed_vector_base<Tp, sizeof(Tp)>
    {
      typedef __packed_vector_base<Tp, sizeof(Tp)> base_type;
      typedef typename Alloc::template rebind<byte_type>::other  byte_alloc_type;
      
      typedef std::vector<byte_type,  byte_alloc_type>  buffer_type;
      typedef utils::succinct_vector<byte_alloc_type>   index_type;
      
      typedef boost::iostreams::filtering_ostream ostream_type;
      
      impl() {}
      ~impl() { close(); }
      
      std::streamsize write(const char_type* s, std::streamsize n);
      void open(const path_type& file);
      bool is_open() const;
      void close();
      
      buffer_type buffer;
      buffer_type buffer_coded;
      size_type   offset;
      size_type   size;
      size_type   size_coded;
      
      boost::shared_ptr<ostream_type> os_data;
      index_type index;
      
      path_type path;
    };
    
    boost::shared_ptr<impl> pimpl;
  };
  
  template <typename Tp,typename Alloc>
  inline
  bool packed_sink<Tp,Alloc>::impl::is_open() const { return os_data; }
  
  template <typename Tp,typename Alloc>
  inline
  std::streamsize packed_sink<Tp,Alloc>::impl::write(const char_type* s, std::streamsize n)
  {
    if (! is_open()) return -1;
    
    const std::streamsize copy_size = std::min(std::streamsize(buffer.size() - offset), n);
    
    std::copy(s, s + copy_size, buffer.begin() + offset);
    offset += copy_size;
    
    if (offset == buffer.size()) {
      const value_type* first = reinterpret_cast<const value_type*>(&(*buffer.begin()));
      const value_type* last = reinterpret_cast<const value_type*>(&(*buffer.end()));
      
      for (const value_type* iter = first; iter != last; ++ iter, ++ size)
	base_type::encode(buffer_coded, size_coded, index, *iter);
      
      if (! buffer_coded.empty() && (sizeof(Tp) > 1 || (size_coded & 0x01) == 0)) {
	os_data->write((char*) &(*buffer_coded.begin()), buffer_coded.size());
	buffer_coded.clear();
	size_coded = 0;
      }
      
      offset = 0;
    }
    
    return copy_size;
  }
  
  template <typename Tp,typename Alloc>
  inline
  void packed_sink<Tp,Alloc>::impl::close()
  {
    typedef utils::repository repository_type;

    if (is_open()) {
      if (offset > 0) {
	if (offset % sizeof(value_type) != 0)
	  throw std::runtime_error("invalid write to this device...");
	
	buffer.resize(offset);
	
	const value_type* first = reinterpret_cast<const value_type*>(&(*buffer.begin()));
	const value_type* last = reinterpret_cast<const value_type*>(&(*buffer.end()));
	
	for (const value_type* iter = first; iter != last; ++ iter, ++ size)
	  base_type::encode(buffer_coded, size_coded, index, *iter);
      }
      
      if (! buffer_coded.empty()) {
	os_data->write((char*) &(*buffer_coded.begin()), buffer_coded.size());
	buffer_coded.clear();
	size_coded = 0;
      }
      
      repository_type rep(path, repository_type::read);
      
      index.build();
      index.write(rep.path("index"));
      
      std::ostringstream stream_size;
      std::ostringstream stream_integral_size;
      stream_size << size;
      stream_integral_size << sizeof(Tp);
      
      rep["size"] = stream_size.str();
      rep["integral-size"] = stream_integral_size.str();
    }
    
    buffer.clear();
    buffer_coded.clear();
    
    offset = 0;
    size = 0;
    size_coded = 0;
    
    os_data.reset();
    index.clear();
    
    path = path_type();
  }
  
  template <typename Tp,typename Alloc>
  inline
  void packed_sink<Tp,Alloc>::impl::open(const path_type& file)
  {
    typedef utils::repository repository_type;
    
    close();
    
    repository_type rep(file, repository_type::write);
    rep["type"] = "packed";
    
    os_data.reset(new boost::iostreams::filtering_ostream());
    os_data->push(boost::iostreams::file_sink(rep.path("data").file_string()), 1024 * 1024);
    
    index.clear();
    
    buffer.reserve(1024 * sizeof(Tp));
    buffer.resize(1024 * sizeof(Tp));
    buffer_coded.clear();
    
    offset = 0;
    size = 0;
    size_coded = 0;
    
    path = file;
  }
};

#endif
