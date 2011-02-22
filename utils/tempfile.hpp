// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// temporry file management
//

#ifndef __UTILS_TEMPFILE__HPP__
#define __UTILS_TEMPFILE__HPP__ 1

#include <sys/types.h>
#include <sys/stat.h>

#include <unistd.h>
#include <errno.h>
#include <cstdlib>
#include <cstring>

#include <string>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <boost/filesystem.hpp>

//
// do we remove singnal installer...?
// we will install signal when tempfile object is created (as static?)
//

namespace utils
{
  class tempfile
  {
  public:
    typedef boost::filesystem::path   path_type;
    
  public:
    static void insert(const path_type& path);
    
    static void erase(const path_type& path);
    
    static path_type tmp_dir();
    
    static path_type file_name(const std::string& file)
    {
      std::vector<char, std::allocator<char> > buffer(file.size() + 1, 0);
      std::copy(file.begin(), file.end(), buffer.begin());
      
      int fd = ::mkstemp(&(*buffer.begin()));
      if (fd < 0)
	throw std::runtime_error(std::string("mkstemp failure: ") + ::strerror(errno));
      ::close(fd);
      
      return path_type(std::string(buffer.begin(), buffer.begin() + file.size()));
    }
    
    static path_type directory_name(const std::string& dir)
    {
      std::vector<char, std::allocator<char> > buffer(dir.size() + 1, 0);
      std::copy(dir.begin(), dir.end(), buffer.begin());
      
      char* tmp = ::mkdtemp(&(*buffer.begin()));
      if (! tmp)
	throw std::runtime_error(std::string("mkdtemp failure: ") + ::strerror(errno));
      
      return path_type(std::string(buffer.begin(), buffer.begin() + dir.size()));
    }
    
#if BOOST_FILESYSTEM_VERSION == 2
    static path_type file_name(const path_type& file) { return file_name(file.file_string()); }
    static path_type directory_name(const path_type& file) { return directory_name(file.file_string()); }
#else
    static path_type file_name(const path_type& file) { return file_name(file.string()); }
    static path_type directory_name(const path_type& file) { return directory_name(file.string()); }
#endif

    static void permission(const path_type& path)
    {
      struct stat buf;
      
#if BOOST_FILESYSTEM_VERSION == 2
      if (::stat(path.file_string().c_str(), &buf) != 0) return;
      
      ::chmod(path.file_string().c_str(), buf.st_mode | S_IRUSR | S_IRGRP | S_IROTH);
#else
      if (::stat(path.string().c_str(), &buf) != 0) return;
      
      ::chmod(path.string().c_str(), buf.st_mode | S_IRUSR | S_IRGRP | S_IROTH);
#endif
    }
  };
};

#endif
