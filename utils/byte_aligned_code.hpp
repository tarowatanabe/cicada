// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__BYTE_ALIGNED_CODE__HPP__
#define __UTILS__BYTE_ALIGNED_CODE__HPP__ 1

#include <stdint.h>
#include <sys/types.h>

namespace utils
{
  template <typename Tp, size_t ByteSize>
  struct __byte_aligned_code {};
  
  template <typename Tp>
  struct __byte_aligned_code<Tp, 1>
  {
    static inline
    size_t byte_size(Tp x) { return 1; }
    
    static inline
    size_t decode(Tp& x, const char* buf)
    {
      x = *buf;
      return 1;
    }

    static inline
    size_t encode(Tp x, char* buf)
    {
      (*buf) = x;
      return 1;
    }
  };
  
  template <typename Tp>
  struct __byte_aligned_code<Tp, 2>
  {
    static inline
    size_t byte_size(Tp x)  { return __byte_size((uint16_t) x); }
    
    static const uint16_t mask1 = ~uint16_t(0x7f);
    static const uint16_t mask2 = uint16_t(mask1 << 7);

    static inline
    size_t __byte_size(uint16_t x)
    {
      return 1 + bool(x & mask1) + bool(x & mask2);
    }
    
    static inline
    size_t decode(Tp& x, const char* buf)
    {
      x = Tp(0);
      size_t __size = 0;
      for (/**/; __size < 3; ++ __size, ++ buf) {
	x |= Tp(*buf & 0x7f) << (__size * 7);
	if (! (*buf & 0x80))
	  return __size + 1;
      }
      return __size;
    }

    static inline
    size_t encode(Tp x, char* buf)
    {
      const size_t __size = byte_size(x);
      switch (__size) {
      case 3: buf[2] = ((x >> 14) & 0x7f) | 0x80;
      case 2: buf[1] = ((x >>  7) & 0x7f) | 0x80;
      case 1: buf[0] = ((x      ) & 0x7f) | 0x80;
      }
      buf[__size - 1] &= 0x7f;
      return __size;
    }
  };
  
  template <typename Tp>
  struct __byte_aligned_code<Tp, 4>
  {
    static inline
    size_t byte_size(Tp x)  { return __byte_size((uint32_t) x); }
    
    static const uint32_t mask1 = ~uint32_t(0x7f);
    static const uint32_t mask2 = mask1 << 7;
    static const uint32_t mask3 = mask2 << 7;
    static const uint32_t mask4 = mask3 << 7;
    
    static inline
    size_t __byte_size(uint32_t x)
    {
      return 1 + bool(x & mask1) + bool(x & mask2) + bool(x & mask3) + bool(x & mask4);
    }
    
    static inline
    size_t decode(Tp& x, const char* buf)
    {
      x = Tp(0);
      size_t __size = 0;
      for (/**/; __size < 5; ++ __size, ++ buf) {
	x |= Tp(*buf & 0x7f) << (__size * 7);
	if (! (*buf & 0x80))
	  return __size + 1;
      }
      return __size;
    }
    
    static inline
    size_t encode(Tp x, char* buf)
    {
      const uint32_t mask = 0x7f;
      const size_t __size = byte_size(x);
      switch (__size) {
      case 5: buf[4] = ((x >> 28) & mask) | 0x80;
      case 4: buf[3] = ((x >> 21) & mask) | 0x80;
      case 3: buf[2] = ((x >> 14) & mask) | 0x80;
      case 2: buf[1] = ((x >>  7) & mask) | 0x80;
      case 1: buf[0] = ((x      ) & mask) | 0x80;
      }
      buf[__size - 1] &= 0x7f;
      return __size;
    }
  };
  
  template <typename Tp>
  struct __byte_aligned_code<Tp, 8>
  {
    static inline
    size_t byte_size(Tp x)  { return __byte_size((uint64_t) x); }
    
    static const uint64_t mask1 = ~uint64_t(0x7f);
    static const uint64_t mask2 = mask1 << 7;
    static const uint64_t mask3 = mask2 << 7;
    static const uint64_t mask4 = mask3 << 7;
    static const uint64_t mask5 = mask4 << 7;
    static const uint64_t mask6 = mask5 << 7;
    static const uint64_t mask7 = mask6 << 7;
    static const uint64_t mask8 = mask7 << 7;
    static const uint64_t mask9 = mask8 << 7;

    static inline
    size_t __byte_size(uint64_t x)
    {
      return (1 
	      + bool(x & mask1) 
	      + bool(x & mask2)
	      + bool(x & mask3)
	      + bool(x & mask4)
	      + bool(x & mask5)
	      + bool(x & mask6)
	      + bool(x & mask7)
	      + bool(x & mask8)
	      + bool(x & mask9));
    }

    static inline
    size_t decode(Tp& x, const char* buf)
    {
      x = Tp(0);
      size_t __size = 0;
      for (/**/; __size < 10; ++ __size, ++ buf) {
	x |= Tp(*buf & 0x7f) << (__size * 7);
	if (! (*buf & 0x80))
	  return __size + 1;
      }
      return __size;
    }
    
    static inline
    size_t encode(Tp x, char* buf)
    {
      const uint64_t uvalue = (uint64_t) x;
      const uint64_t mask = 0x7f;
      const size_t __size = byte_size(uvalue);
      switch (__size) {
      case 10: buf[9] = ((uvalue >> 63) & mask) | 0x80;
      case  9: buf[8] = ((uvalue >> 56) & mask) | 0x80;
      case  8: buf[7] = ((uvalue >> 49) & mask) | 0x80;
      case  7: buf[6] = ((uvalue >> 42) & mask) | 0x80;
      case  6: buf[5] = ((uvalue >> 35) & mask) | 0x80;
      case  5: buf[4] = ((uvalue >> 28) & mask) | 0x80;
      case  4: buf[3] = ((uvalue >> 21) & mask) | 0x80;
      case  3: buf[2] = ((uvalue >> 14) & mask) | 0x80;
      case  2: buf[1] = ((uvalue >>  7) & mask) | 0x80;
      case  1: buf[0] = ((uvalue      ) & mask) | 0x80;
      }
      buf[__size - 1] &= 0x7f;
      return __size;
    }
  };
  
  template <typename Tp>
  size_t byte_aligned_encode(const Tp& x, char* buf)
  {
    return __byte_aligned_code<Tp, sizeof(Tp)>::encode(x, buf);
  }
  
  template <typename Tp>
  size_t byte_aligned_decode(Tp& x, const char* buf)
  {
    return __byte_aligned_code<Tp, sizeof(Tp)>::decode(x, buf);
  }
};

#endif
