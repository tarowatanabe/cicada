// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__COMPRESS_STREAM__HPP__
#define __UTILS__COMPRESS_STREAM__HPP__ 1

#include <cstring>
#include <iostream>
#include <fstream>

#include <unistd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>

#include <boost/version.hpp>
#include <boost/filesystem.hpp>

namespace utils
{
  namespace impl
  {
    
    typedef enum {
      COMPRESS_STREAM_GZIP,
      COMPRESS_STREAM_BZIP,
      COMPRESS_STREAM_UNKNOWN
    } compress_format_type;
    
    inline compress_format_type compress_iformat (const std::string& filename)
    {
      char buffer[8];
      
      ::memset(buffer, 0, sizeof(char)*3);
      std::ifstream ifs(filename.c_str());
      ifs.read((char*) buffer, sizeof(char)*3);
      
      if (buffer[0] == '\037' && buffer[1] == '\213')
	return COMPRESS_STREAM_GZIP;
      else if (buffer[0] == 'B' && buffer[1] == 'Z' && buffer[2] == 'h')
	return COMPRESS_STREAM_BZIP;
      else
	return COMPRESS_STREAM_UNKNOWN;
    }
    inline compress_format_type compress_iformat(const boost::filesystem::path& path) { return compress_iformat(path.file_string()); }
    
    inline compress_format_type compress_oformat (const std::string& filename)
    {
      if (filename.size() > 3 && strncmp(&(filename.c_str()[filename.size() - 3]), ".gz", 3) == 0)
	return COMPRESS_STREAM_GZIP;
      else if (filename.size() > 4 && strncmp(&(filename.c_str()[filename.size() - 4]), ".bz2", 4) == 0)
	return COMPRESS_STREAM_BZIP;
      else
	return COMPRESS_STREAM_UNKNOWN;
    }
    inline compress_format_type compress_oformat(const boost::filesystem::path& path) { return compress_oformat(path.file_string()); }    
  };


  template <typename Stream>
  inline
  Stream& push_compress_ostream(Stream& os,
				const boost::filesystem::path& path,
				size_t buffer_size = 4096)
  {
    if (path.file_string() == "-")  {
#if BOOST_VERSION >= 104400
      os.push(boost::iostreams::file_descriptor_sink(::dup(STDOUT_FILENO), boost::iostreams::close_handle), buffer_size);
#else
      os.push(boost::iostreams::file_descriptor_sink(::dup(STDOUT_FILENO), true), buffer_size);
#endif
    } else {
      switch (impl::compress_oformat(path)) {
      case impl::COMPRESS_STREAM_GZIP:
	os.push(boost::iostreams::gzip_compressor());
	break;
      case impl::COMPRESS_STREAM_BZIP:
	os.push(boost::iostreams::bzip2_compressor());
	break;
      default: break;
      }
      os.push(boost::iostreams::file_sink(path.file_string(), std::ios_base::out | std::ios_base::trunc), buffer_size);
    }
    return os;
  }
  
  template <typename Stream>
  inline
  Stream& push_compress_istream(Stream& is,
				const boost::filesystem::path& path,
				size_t buffer_size = 4096)
  {
    if (path.file_string() == "-")  {
#if BOOST_VERSION >= 104400
      is.push(boost::iostreams::file_descriptor_source(::dup(STDIN_FILENO), boost::iostreams::close_handle), buffer_size);
#else
      is.push(boost::iostreams::file_descriptor_source(::dup(STDIN_FILENO), true), buffer_size);
#endif
    } else {
      if (boost::filesystem::is_regular_file(path)) {
	switch (impl::compress_iformat(path)) {
	case impl::COMPRESS_STREAM_GZIP:
	  is.push(boost::iostreams::gzip_decompressor());
	  break;
	case impl::COMPRESS_STREAM_BZIP:
	  is.push(boost::iostreams::bzip2_decompressor());
	  break;
	default: break;
	}
      }
      is.push(boost::iostreams::file_source(path.file_string()), buffer_size);
    }
    return is;
  }
  
  class compress_ostream : public boost::iostreams::filtering_ostream
  {
  public:
    typedef boost::filesystem::path path_type;
  private:
    typedef boost::iostreams::filtering_ostream stream_type;
    
  public:
    compress_ostream(const path_type& path,
		     size_t buffer_size = 4096)
    {
      push_compress_ostream(__stream(), path, buffer_size);
    }
    
  private:
    stream_type& __stream() { return static_cast<stream_type&>(*this); }
  };
  
  class compress_istream : public boost::iostreams::filtering_istream
  {
  public:
    typedef boost::filesystem::path path_type;

  private:
    typedef boost::iostreams::filtering_istream stream_type;
    
  public:
    compress_istream(const path_type& path,
		     size_t buffer_size = 4096)
    {
      push_compress_istream(__stream(), path, buffer_size);
    }
    
  private:
    stream_type& __stream() { return static_cast<stream_type&>(*this); }
  };
  
};

#endif
