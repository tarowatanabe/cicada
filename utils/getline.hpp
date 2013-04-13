//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//
// -*- mode: c++ -*-

#ifndef __UTILS__GETLINE__HPP__
#define __UTILS__GETLINE__HPP__ 1

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace utils
{
  inline
  std::istream& getline(std::istream& is, std::string& line)
  {
    char buf[4096];
    std::vector<char, std::allocator<char> > buffer;

    line.clear();
    
    do {
      is.clear();
      is.getline(buf, 4096);
      buffer.insert(buffer.end(), buf, buf + std::strlen(buf));
    } while (! is.eof() && is.fail());
    
    if (! is.eof())
      is.clear();
    
    line.assign(buffer.begin(), buffer.end());
    
    return is;
  }
  
};

#endif
