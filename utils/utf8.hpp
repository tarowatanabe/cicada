// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__UTF8__HPP__
#define __UTILS__UTF8__HPP__ 1

#include<stdexcept>

namespace utils
{
  inline
  bool utf8_start(const unsigned char c)
  {
    return ((c & 0x80) == 0) || ((c & 0xC0) == 0xC0);
  }
  
  inline
  bool utf8_other(const unsigned char c)
  {
    return (c & 0xc0) == 0x80;
  }
  
  inline
  size_t utf8_size(const unsigned char c)
  {
#if 0
    if ((c & 0x80) == 0) 
      return 1;
    else if ((c & 0xE0) == 0xC0)
      return 2;
    else if ((c & 0xF0) == 0xE0)
      return 3;
    else if ((c & 0xF8) == 0xF0)
      return 4;
    else if ((c & 0xFC) == 0xF8)
      return 5;
    else if ((c & 0xFE) == 0xFC)
      return 6;
    else
      throw std::runtime_error("invlaid utf-8 sequence");
#endif

#if 1
    static const char bytes_utf8[256] = {
      1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,
      2,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,
      3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3, 4,4,4,4,4,4,4,4, 5,5,5,5,6,6,-1,-1};
    
    return bytes_utf8[c];
#endif
  }
};

#endif
