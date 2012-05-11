// -*- mode: c++ -*-

#ifndef __CICADA_OUTPUT_IMPL__HPP__
#define __CICADA_OUTPUT_IMPL__HPP__ 1

#include <unistd.h>

#include <boost/filesystem.hpp>

#include "utils/filesystem.hpp"

inline
void prepare_directory(const boost::filesystem::path& path)
{
  if (boost::filesystem::exists(path) && ! boost::filesystem::is_directory(path))
    utils::filesystem::remove_all(path);
  
  if (! boost::filesystem::exists(path))
    boost::filesystem::create_directories(path);
  
  boost::filesystem::directory_iterator iter_end;
  for (boost::filesystem::directory_iterator iter(path); iter != iter_end; ++ iter)
    utils::filesystem::remove_all(*iter);  
  
  ::sync();
}

#endif
