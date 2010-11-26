// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__FILESYSTEM__HPP__
#define __UTILS__FILESYSTEM__HPP__ 1


#include <dirent.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cstdio>
#include <vector>

#include <boost/filesystem.hpp>

// bugfixed version of remove-all, buffer optimized copy_file and copy_files...

namespace utils
{
  namespace filesystem
  {

    inline
    void copy_file(const boost::filesystem::path& path_from,
                   const boost::filesystem::path& path_to)
    {
      // 4M bytes buffer...
      std::vector<char> buffer(4 * 1024 * 1024);

      int infile = 0;
      int outfile = 0;
      struct stat stat_from;

      const std::string file_from = path_from.native_file_string();
      const std::string file_to = path_to.native_file_string();
      
      if (:: stat(file_from.c_str(), &stat_from) != 0
          || (infile = ::open(file_from.c_str(), O_RDONLY)) < 0
          || (outfile = ::open(file_to.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_EXCL, stat_from.st_mode)) < 0) {

        if (infile >= 0) ::close(infile);
        throw std::runtime_error("copy file failed");
      }
      
      ssize_t sz = 0;
      ssize_t sz_read = 1;
      ssize_t sz_write = 0;
      
      while (sz_read > 0) {
        do
          sz_read = ::read(infile, &(*buffer.begin()), buffer.size());
        while (sz_read < 0 && errno == EINTR);

        if (sz_read <= 0) break;

        sz_write = 0;
        do {
          if ((sz = ::write(outfile, &(*(buffer.begin() + sz_write)), sz_read - sz_write)) < 0) {
            if (errno == EINTR)
              continue;
            else {
              sz_read = sz;
              break;
            }
          }
          sz_write += sz;
        } while (sz_write < sz_read);
      }
      if ( ::close( infile) < 0 ) sz_read = -1;
      if ( ::close( outfile) < 0 ) sz_read = -1;

      if ( sz_read < 0 )
        throw std::runtime_error("copy file failed");
    }

    inline
    void copy_files(const boost::filesystem::path& path_from,
		    const boost::filesystem::path& path_to)
    {
      if (path_from.empty() || path_to.empty())
	throw std::runtime_error("empty path");
      
      if (! boost::filesystem::exists(path_from) || ! boost::filesystem::exists(path_to))
	throw std::runtime_error("no path");
      
      if (! boost::filesystem::is_directory(path_to))
	throw std::runtime_error("destination is not a directory");
      
      boost::filesystem::path destination = path_to / path_from.leaf();
      
      if (! boost::filesystem::is_directory(path_from))
	boost::filesystem::copy_file(path_from, destination);
      else {
	boost::filesystem::create_directory(destination);
	
	boost::filesystem::directory_iterator iter_end;
	for (boost::filesystem::directory_iterator iter(path_from); iter != iter_end; ++ iter)
	  copy_files(*iter, destination);
      }
    }
    
    inline
    bool __exists(const boost::filesystem::path& path)
    {
      if (path.empty()) return false;
      struct stat statbuf;
      if (::lstat(path.file_string().c_str(), &statbuf) >= 0)
	return true;
      else {
	errno = 0; return false;
      }
    }
    
    inline
    bool __is_directory(const boost::filesystem::path& path)
    {
      if (path.empty()) return false;
      struct stat statbuf;
      if (::lstat(path.file_string().c_str(), &statbuf) >= 0)
	return ! S_ISLNK(statbuf.st_mode) && S_ISDIR(statbuf.st_mode);
      else {
	errno = 0; return false;
      }
    }
    
    inline
    void remove_all(const boost::filesystem::path& path)
    {
      if (__is_directory(path)) {
	std::vector<boost::filesystem::path> paths;
	
	DIR* dirp = (path.empty() ? ::opendir(".") : ::opendir(path.directory_string().c_str()));
	if (dirp) {
	  struct dirent* dp;
	  do {
	    errno = 0;
	    if ((dp = readdir(dirp)) != 0) {
	      const std::string path_leaf(dp->d_name);
	      if (path_leaf != "." && path_leaf != "..")
		paths.push_back(path / path_leaf);
	    }
	  } while (dp);
	  ::closedir(dirp);
	}
	std::for_each(paths.begin(), paths.end(), utils::filesystem::remove_all);
      }
      if (__exists(path))
	::remove(path.file_string().c_str());
      errno = 0;
    }
  };
};

#endif
