// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

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
    
    line.clear();
    do {
      is.clear();
      is.getline(buf, 4096);
      line.append(buf);
    } while (! is.eof() && is.fail());
    
    if (! is.eof())
      is.clear();
    
    if (line.capacity() > (1 << 30))
      std::string(line).swap(line);
    
    return is;
  }
  
};

#endif
