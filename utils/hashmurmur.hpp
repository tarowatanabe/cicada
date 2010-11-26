// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__HASHMURMUR__HPP__
#define __UTILS__HASHMURMUR__HPP__ 1

#include <cstddef>
#include <stdint.h>

// murmurhash2 by Austin Appleby
// we assume only POD data... but potentially optimized with template metaprogramming
// if we know the size in advance...

namespace utils
{
  
  template <typename Tp, size_t _Remain>
  struct __static_hashmurmur_remain {};
  
  template <typename Tp>
  struct __static_hashmurmur_remain<Tp, 3>
  {
    static inline
    void remain(const uint8_t* data, Tp& h, const Tp m)
    {
      h ^= data[2] << 16;
      h ^= data[1] << 8;
      h ^= data[0];
      h *= m;
    }
  };
  
  template <typename Tp>
  struct __static_hashmurmur_remain<Tp, 2>
  {
    static inline
    void remain(const uint8_t* data, Tp& h, const Tp m)
    {
      h ^= data[1] << 8;
      h ^= data[0];
      h *= m;
    }
  };
  
  template <typename Tp>
  struct __static_hashmurmur_remain<Tp, 1>
  {
    static inline
    void remain(const uint8_t* data, Tp& h, const Tp m)
    {
      h ^= data[0];
      h *= m;
    }
  };
  
  template <typename Tp>
  struct __static_hashmurmur_remain<Tp, 0>
  {
    static inline
    void remain(const uint8_t* data, Tp& h, const Tp m) {}
  };
    
  
  
  template <size_t _Length, size_t _Remain, size_t _Size>
  struct __static_hashmurmur_impl {};
  
  template <size_t _Length, size_t _Remain>
  struct __static_hashmurmur_impl<_Length, _Remain, 4>
  {    
    static inline
    uint32_t hash(const uint8_t* data, uint32_t seed)
    {
      const uint32_t m = 0x5bd1e995;
      const int32_t r  = 24;
      
      size_t length = _Length;
      uint32_t h = seed ^ _Length;
      
      for (/**/; length >= 4; length -= 4, data += 4) {
	uint32_t k = *((uint32_t*) data);
	
	k *= m; 
	k ^= k >> r; 
	k *= m;
	  
	h *= m; 
	h ^= k;
      }
      
      // Handle the last few bytes of the input array
      __static_hashmurmur_remain<uint32_t, _Remain>::remain(data, h, m);
      
      // Do a few final mixes of the hash to ensure the last few
      // bytes are well-incorporated.
      h ^= h >> 13;
      h *= m;
      h ^= h >> 15;
	
      return h;
    }
  };
  
  template <size_t _Length>
  struct __static_hashmurmur_impl<_Length, 0, 4>
  {
    static inline
    uint32_t hash(const uint8_t* __data, uint32_t seed)
    {
      const uint32_t m   = 0x5bd1e995;
      const int32_t r = 24;
      
      const uint32_t* data = ((const uint32_t*) __data);
      uint32_t h = seed ^ _Length;
      
      for (size_t i = 0; i < (_Length >> 2); ++ i, ++ data) {
	uint32_t k = *data;
	  
	k *= m; 
	k ^= k >> r; 
	k *= m; 
	  
	h *= m; 
	h ^= k;
      }
	
      // Do a few final mixes of the hash to ensure the last few
      // bytes are well-incorporated.
      h ^= h >> 13;
      h *= m;
      h ^= h >> 15;
	
      return h;
    }
  };
  
  template <size_t _Length, size_t _Remain>
  struct __static_hashmurmur_impl<_Length, _Remain, 8>
  {    
    static inline
    uint64_t hash(const uint8_t* data, uint64_t seed)
    {
      const uint32_t m   = 0x5bd1e995;
      const uint64_t m64 = 0x5bd1e995;
      const int32_t r = 24;
      
      size_t length = _Length;
      uint64_t h = seed ^ _Length;
	
      for (/**/; length >= 4; length -= 4, data += 4) {
	uint32_t k = *((uint32_t*) data);
	  
	k *= m; 
	k ^= k >> r; 
	k *= m;
	
	h *= m64; 
	h ^= k;
      }
      
      // Handle the last few bytes of the input array
      __static_hashmurmur_remain<uint64_t, _Remain>::remain(data, h, m64);
      
      // Do a few final mixes of the hash to ensure the last few
      // bytes are well-incorporated.
      h ^= h >> 13;
      h *= m64;
      h ^= h >> 15;
      
      return h;
    }
  };
  
  template <size_t _Length>
  struct __static_hashmurmur_impl<_Length, 0, 8>
  {
    static inline
    uint64_t hash(const uint8_t* __data, uint64_t seed)
    {
      const uint32_t m   = 0x5bd1e995;
      const uint64_t m64 = 0x5bd1e995;
      const int32_t r = 24;
	
      const uint32_t* data = ((const uint32_t*) __data);
      uint64_t h = seed ^ _Length;
	
      for (size_t i = 0; i < (_Length >> 2); ++ i, ++ data) {
	uint32_t k = *data;
	  
	k *= m; 
	k ^= k >> r; 
	k *= m; 
	  
	h *= m64; 
	h ^= k;
      }
	
      // Do a few final mixes of the hash to ensure the last few
      // bytes are well-incorporated.
      h ^= h >> 13;
      h *= m64;
      h ^= h >> 15;
	
      return h;
    }
  };
  
  template <size_t _Length, size_t _Size>
  struct __static_hashmurmur {};
  
  template <size_t _Length>
  struct __static_hashmurmur<_Length, 4>
  {
    static inline
    uint32_t hash(const uint8_t* data, uint32_t seed)
    {
      return __static_hashmurmur_impl<_Length, _Length & 0x03, 4>::hash(data, seed);
    }
  };
  
  template <size_t _Length>
  struct __static_hashmurmur<_Length, 8>
  {
    static inline
    uint32_t hash(const uint8_t* data, uint32_t seed)
    {
      return __static_hashmurmur_impl<_Length, _Length & 0x03, 8>::hash(data, seed);
    }
  };
  
  template <size_t _Sizse>
  struct __dynamic_hashmurmur {};
  
  template <>
  struct __dynamic_hashmurmur<4>
  {
    template <typename Iterator>
    static inline
    uint32_t hash(Iterator first, Iterator last, uint32_t seed)
    {
      const uint32_t m = 0x5bd1e995;
      const int32_t r = 24;
      
      uint32_t h = seed ^ (last - first);
      size_t length = last - first;
      
      while (length >= 4) {
	uint32_t k = *((uint32_t*) first);
	
	k *= m;
	k ^= k >> r; 
	k *= m;
	
	h *= m; 
	h ^= k;
	
	first += 4;
	length -= 4;
      }
      
      // Handle the last few bytes of the input array
      switch (length) {
      case 3: h ^= *(first + 2) << 16;
      case 2: h ^= *(first + 1) << 8;
      case 1: h ^= *(first + 0);
	h *= m;
      };
      
      // Do a few final mixes of the hash to ensure the last few
      // bytes are well-incorporated.
      h ^= h >> 13;
      h *= m;
      h ^= h >> 15;
      
      return h;
    }
  };
  
  template <>
  struct __dynamic_hashmurmur<8>
  {
    template <typename Iterator>
    static inline
    uint64_t hash(Iterator first, Iterator last, uint64_t seed)
    {
      const uint32_t m   = 0x5bd1e995;
      const uint64_t m64 = 0x5bd1e995;
      const int32_t r = 24;
      
      uint64_t h = seed ^ (last - first);
      size_t length = last - first;
      
      while (length >= 4) {
	uint32_t k = *((uint32_t*) first);
	
	k *= m; 
	k ^= k >> r; 
	k *= m; 
	
	h *= m64; 
	h ^= k;
	
	first += 4;
	length -= 4;
      }
      
      // Handle the last few bytes of the input array
      switch (length) {
      case 3: h ^= *(first + 2) << 16;
      case 2: h ^= *(first + 1) << 8;
      case 1: h ^= *(first + 0);
	h *= m64;
      };
      
      // Do a few final mixes of the hash to ensure the last few
      // bytes are well-incorporated.
      h ^= h >> 13;
      h *= m64;
      h ^= h >> 15;
      
      return h;
    }
  };
  
  template <typename _Value>
  struct hashmurmur
  {
    
  public:
    template <typename Tp, size_t N>
    _Value operator()(const Tp (&x)[N], _Value seed=0) const
    { 
      return __static_hashmurmur<sizeof(Tp) * N, sizeof(_Value)>::hash((const uint8_t*) x, seed);
    }
    
    template <typename Tp>
    _Value operator()(const Tp& x, _Value seed=0) const
    {
      return __static_hashmurmur<sizeof(Tp), sizeof(_Value)>::hash((const uint8_t*) &x, seed);
    }
    
    template <typename Iterator>
    _Value operator()(Iterator first, Iterator last, _Value seed) const
    {
      return __dynamic_hashmurmur<sizeof(_Value)>::hash((const uint8_t*) &(*first), (const uint8_t*) &(*last), seed);
    }
  };
  
};

#endif
