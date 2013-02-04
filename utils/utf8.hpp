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
  
  template <typename Buf>
  inline
  int utf8_code(Buf buf)
  {
    typedef unsigned int code_type;
    
    switch (utf8_size(*buf)) {
    case 1: return *buf;
    case 2: return ((code_type(*buf & 0x1F) << 6)
		    | (code_type(*(buf + 1) & 0x3F)));
    case 3: return ((code_type(*buf & 0x0F) << 12)
		    | (code_type(*(buf + 1) & 0x3F) <<  6)
		    | (code_type(*(buf + 2) & 0x3F)));
    case 4: return ((code_type(*buf & 0x07) << 18)
		    | (code_type(*(buf + 1) & 0x3F) << 12)
		    | (code_type(*(buf + 2) & 0x3F) << 6)
		    | (code_type(*(buf + 3) & 0x3F)));
    case 5: return ((code_type(*buf & 0x03) << 24)
		    | (code_type(*(buf + 1) & 0x3F) << 18)
		    | (code_type(*(buf + 2) & 0x3F) << 12)
		    | (code_type(*(buf + 3) & 0x3F) << 6)
		    | (code_type(*(buf + 4) & 0x3F)));
    case 6: return ((code_type(*buf & 0x01) << 30)
		    | (code_type(*(buf + 1) & 0x3F) << 24)
		    | (code_type(*(buf + 2) & 0x3F) << 18)
		    | (code_type(*(buf + 3) & 0x3F) << 12)
		    | (code_type(*(buf + 4) & 0x3F) << 6)
		    | (code_type(*(buf + 5) & 0x3F)));
    default:
      throw std::runtime_error("invlaid utf-8");
    }
  }

  template <typename Sink>
  inline
  Sink& utf8_sink(Sink& sink, unsigned int code)
  {
    if (code <= 0x7f) {
      sink.put(code);
    } else if (code <= 0x7ff) {
      sink.put(((code >> 6) & 0x1f) | 0xc0);
      sink.put(((code     ) & 0x3f) | 0x80);
    } else if (code <= 0xffff) {
      sink.put(((code >> 12) & 0x0f) | 0xe0);
      sink.put(((code >>  6) & 0x3f) | 0x80);
      sink.put(((code      ) & 0x3f) | 0x80);
    } else if (code <= 0x1fffff) {
      sink.put(((code >> 18) & 0x07) | 0xf0);
      sink.put(((code >> 12) & 0x3f) | 0x80);
      sink.put(((code >>  6) & 0x3f) | 0x80);
      sink.put(((code      ) & 0x3f) | 0x80);
    } else if (code <= 0x3ffffff) {
      sink.put(((code >> 24) & 0x03) | 0xf8);
      sink.put(((code >> 18) & 0x3f) | 0x80);
      sink.put(((code >> 12) & 0x3f) | 0x80);
      sink.put(((code >>  6) & 0x3f) | 0x80);
      sink.put(((code      ) & 0x3f) | 0x80);
    } else if (code <= 0x7fffffff) {
      sink.put(((code >> 30) & 0x01) | 0xfc);
      sink.put(((code >> 24) & 0x3f) | 0x80);
      sink.put(((code >> 18) & 0x3f) | 0x80);
      sink.put(((code >> 12) & 0x3f) | 0x80);
      sink.put(((code >>  6) & 0x3f) | 0x80);
      sink.put(((code      ) & 0x3f) | 0x80);
    } else
      throw std::runtime_error("out of range");
    
    return sink;
  }
};

#endif
