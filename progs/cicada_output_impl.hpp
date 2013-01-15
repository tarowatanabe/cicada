// -*- mode: c++ -*-

#ifndef __CICADA_OUTPUT_IMPL__HPP__
#define __CICADA_OUTPUT_IMPL__HPP__ 1

#include <unistd.h>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

#include "utils/filesystem.hpp"

inline
void prepare_directory(const boost::filesystem::path& path)
{
  if (! boost::filesystem::exists(path))
    boost::filesystem::create_directories(path);
  else if (! boost::filesystem::is_directory(path)) {
    utils::filesystem::remove_all(path);
    boost::filesystem::create_directories(path);
  }

  while (! boost::filesystem::exists(path)) {
    ::sync();
    boost::thread::yield();
  }
  
  boost::filesystem::directory_iterator iter_end;
  for (boost::filesystem::directory_iterator iter(path); iter != iter_end; ++ iter) {
    const boost::filesystem::path subdir = *iter;
    
    utils::filesystem::remove_all(subdir);
    
    while (boost::filesystem::exists(subdir)) {
      ::sync();
      boost::thread::yield();
    }
  }

  ::sync();
}

#endif
