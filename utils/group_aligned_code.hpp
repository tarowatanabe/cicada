// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__GROUP_ALIGNED_CODE__HPP__
#define __UTILS__GROUP_ALIGNED_CODE__HPP__ 1

#include <stdint.h>
#include <sys/types.h>

namespace utils
{
  
  template <typename Tp, size_t ByteSize>
  struct __group_aligned_code {};
  
  template <typename Tp>
  struct __group_aligned_code<Tp, 4>
  {
    typedef Tp       value_type;
    typedef uint32_t uvalue_type;
    

    static inline
    uint8_t* offsets()
    {
      static uint8_t __table[256 * 5] = {
#include <utils/group_aligned_code_offsets.hpp>
      };
      return __table;
    }
    
    static inline 
    size_t byte_size(const Tp& x)
    {
      return __byte_size((uvalue_type) x);
    }
    
    static inline 
    size_t __byte_size(uvalue_type x)
    {
      return 1 + bool(x & 0xffffff00) + bool(x & 0xffff0000) + bool(x & 0xff000000);
    }
    
    static inline
    size_t decode(Tp& x, const char* buf, size_t pos)
    {
      pos &= 0x03;
      
      const size_t first = 1 + offsets()[uint8_t(*buf) * 5 + pos];
      const size_t last  = 1 + offsets()[uint8_t(*buf) * 5 + pos + 1];

      switch (last - first) {
      case 1: x = Tp(buf[first]) & 0xff; break;
      case 2: x = (((Tp(buf[first]) & 0xff) << 8)
		   | ((Tp(buf[first + 1]) & 0xff))); break;
      case 3: x = (((Tp(buf[first]) & 0xff) << 16)
		   | ((Tp(buf[first + 1]) & 0xff) << 8)
		   | ((Tp(buf[first + 2]) & 0xff))); break;
      case 4: x = (((Tp(buf[first]) & 0xff) << 24)
		   | ((Tp(buf[first + 1]) & 0xff) << 16)
		   | ((Tp(buf[first + 2]) & 0xff) << 8)
		   | ((Tp(buf[first + 3]) & 0xff))); break;
      }
      return last;
    }
    
    static inline
    size_t encode(const Tp& x, char* buf, size_t pos)
    {
      pos &= 0x03;
      
      const size_t __size = byte_size(x);
      
      const size_t first = 1 + offsets()[uint8_t(*buf) * 5 + pos];
      const size_t last = first + __size;
      
      // set buf indicating size of x
      const size_t shift_amount = (3 - pos) * 2;
      
      *buf ^= (*buf ^ ((__size - 1) << shift_amount)) & (0x03 << shift_amount);
      
      switch (__size) {
      case 1: buf[first] = x; break;
      case 2: buf[first] = (x >>  8); buf[first + 1] = x; break;
      case 3: buf[first] = (x >> 16); buf[first + 1] = (x >>  8); buf[first + 2] = x; break;
      case 4: buf[first] = (x >> 24); buf[first + 1] = (x >> 16); buf[first + 2] = (x >> 8); buf[first + 3] = x; break;
      }
      return last;
    }
  };

  template <typename Tp>
  size_t group_aligned_encode(const Tp& x, char* buf, size_t pos)
  {
    return __group_aligned_code<Tp, sizeof(Tp)>::encode(x, buf, pos);
  }
  
  template <typename Tp>
  size_t group_aligned_decode(Tp& x, const char* buf, size_t pos)
  {
    return __group_aligned_code<Tp, sizeof(Tp)>::decode(x, buf, pos);
  }  
};

#endif
