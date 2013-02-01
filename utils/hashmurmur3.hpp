// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//-----------------------------------------------------------------------------
// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.
//
// MIT license
//

#ifndef __UTILS__HASHMURMUR3__HPP__
#define __UTILS__HASHMURMUR3__HPP__ 1

#include <cstddef>
#include <stdint.h>

namespace utils
{

  template <size_t _Size>
  struct __dynamic_hashmurmur3 {};
  
  template <>
  struct __dynamic_hashmurmur3<4>
  {
    static inline
    uint32_t getblock ( const uint32_t * p, int i )
    {
      return p[i];
    }
    
    static inline
    uint32_t rotl32 ( uint32_t x, int8_t r )
    {
      return (x << r) | (x >> (32 - r));
    }

    static inline
    uint32_t fmix ( uint32_t h )
    {
      h ^= h >> 16;
      h *= 0x85ebca6b;
      h ^= h >> 13;
      h *= 0xc2b2ae35;
      h ^= h >> 16;
      
      return h;
    }

    template <typename Iterator>
    static inline
    uint32_t hash(Iterator first, Iterator last, uint32_t seed)
    {
      const uint8_t * data = (const uint8_t*) first;
      size_t len = last - first;
      const ptrdiff_t nblocks = len >> 2;
      
      uint32_t h1 = seed;
      
      const uint32_t c1 = 0xcc9e2d51;
      const uint32_t c2 = 0x1b873593;
      
      //----------
      // body
      
      const uint32_t * blocks = (const uint32_t *)(data + (nblocks << 2));
      
      for (ptrdiff_t i = - nblocks; i; i++) {
	uint32_t k1 = getblock(blocks,i);
	
	k1 *= c1;
	k1 = rotl32(k1,15);
	k1 *= c2;
	
	h1 ^= k1;
	h1 = rotl32(h1,13); 
	h1 = h1*5+0xe6546b64;
      }
      
      //----------
      // tail
      
      const uint8_t * tail = (const uint8_t*)(data + (nblocks << 2));
      
      uint32_t k1 = 0;
      
      switch(len & 3) {
      case 3: k1 ^= tail[2] << 16;
      case 2: k1 ^= tail[1] << 8;
      case 1: k1 ^= tail[0];
	k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1;
      };
      
      //----------
      // finalization
      
      h1 ^= static_cast<uint32_t>(len);
      
      return fmix(h1);
    }
  };
  
  template <>
  struct __dynamic_hashmurmur3<8>
  {
    static inline
    uint64_t getblock ( const uint64_t * p, int i )
    {
      return p[i];
    }
    
    static inline
    uint64_t rotl64 ( uint64_t x, int8_t r )
    {
      return (x << r) | (x >> (64 - r));
    }

    static inline
    uint64_t fmix ( uint64_t k )
    {
      k ^= k >> 33;
      k *= 0xff51afd7ed558ccdULL;
      k ^= k >> 33;
      k *= 0xc4ceb9fe1a85ec53ULL;
      k ^= k >> 33;
      
      return k;
    }
    
    template <typename Iterator>
    static inline
    uint64_t hash(Iterator first, Iterator last, uint64_t seed)
    {
      const uint8_t * data = (const uint8_t*)first;
      size_t len = last - first;
      const size_t nblocks = len >> 4;
      
      uint64_t h1 = seed;
      uint64_t h2 = seed;
      
      const uint64_t c1 = 0x87c37b91114253d5LLU;
      const uint64_t c2 = 0x4cf5ad432745937fLLU;
      
      //----------
      // body
      
      const uint64_t * blocks = (const uint64_t *)(data);
      
      for (size_t i = 0; i != nblocks; i++) {
	uint64_t k1 = getblock(blocks,i*2+0);
	uint64_t k2 = getblock(blocks,i*2+1);
	
	k1 *= c1; k1  = rotl64(k1,31); k1 *= c2; h1 ^= k1;
	
	h1 = rotl64(h1,27); h1 += h2; h1 = h1*5+0x52dce729;
	
	k2 *= c2; k2  = rotl64(k2,33); k2 *= c1; h2 ^= k2;
	
	h2 = rotl64(h2,31); h2 += h1; h2 = h2*5+0x38495ab5;
      }
      
      //----------
      // tail
      
      const uint8_t * tail = (const uint8_t*)(data + (nblocks << 4));
      
      uint64_t k1 = 0;
      uint64_t k2 = 0;
      
      switch(len & 15) {
      case 15: k2 ^= uint64_t(tail[14]) << 48;
      case 14: k2 ^= uint64_t(tail[13]) << 40;
      case 13: k2 ^= uint64_t(tail[12]) << 32;
      case 12: k2 ^= uint64_t(tail[11]) << 24;
      case 11: k2 ^= uint64_t(tail[10]) << 16;
      case 10: k2 ^= uint64_t(tail[ 9]) << 8;
      case  9: k2 ^= uint64_t(tail[ 8]) << 0;
	k2 *= c2; k2  = rotl64(k2,33); k2 *= c1; h2 ^= k2;
	  
      case  8: k1 ^= uint64_t(tail[ 7]) << 56;
      case  7: k1 ^= uint64_t(tail[ 6]) << 48;
      case  6: k1 ^= uint64_t(tail[ 5]) << 40;
      case  5: k1 ^= uint64_t(tail[ 4]) << 32;
      case  4: k1 ^= uint64_t(tail[ 3]) << 24;
      case  3: k1 ^= uint64_t(tail[ 2]) << 16;
      case  2: k1 ^= uint64_t(tail[ 1]) << 8;
      case  1: k1 ^= uint64_t(tail[ 0]) << 0;
	k1 *= c1; k1  = rotl64(k1,31); k1 *= c2; h1 ^= k1;
      };
      
      //----------
      // finalization
      
      h1 ^= len; h2 ^= len;
      
      h1 += h2;
      h2 += h1;
      
      h1 = fmix(h1);
      h2 = fmix(h2);
      
      h1 += h2;
      h2 += h1;
      
      //((uint64_t*)out)[0] = h1; // low
      //((uint64_t*)out)[1] = h2; // high
      
      // Hash128to64 from city hash...
      const uint64_t kMul = 0x9ddfea08eb382d69ULL;
      uint64_t a = (h1 ^ h2) * kMul;
      a ^= (a >> 47);
      uint64_t b = (h2 ^ a) * kMul;
      b ^= (b >> 47);
      b *= kMul;
      return b;
    }
  };

  template <typename _Value>
  struct hashmurmur3
  {
  public:
    template <typename Tp, size_t N>
    _Value operator()(const Tp (&x)[N], _Value seed=0) const
    { 
      //return __static_hashmurmur3<sizeof(Tp) * N, sizeof(_Value)>::hash((const uint8_t*) x, seed);
      
      return __dynamic_hashmurmur3<sizeof(_Value)>::hash((const uint8_t*) &x, ((const uint8_t*) &x) + sizeof(Tp) * N, seed);
    }
    
    template <typename Tp>
    _Value operator()(const Tp& x, _Value seed=0) const
    {
      //return __static_hashmurmur3<sizeof(Tp), sizeof(_Value)>::hash((const uint8_t*) &x, seed);
      
      return __dynamic_hashmurmur3<sizeof(_Value)>::hash((const uint8_t*) &x, ((const uint8_t*) &x) + sizeof(Tp), seed);
    }
    
    template <typename Iterator>
    _Value operator()(Iterator first, Iterator last, _Value seed) const
    {
      return __dynamic_hashmurmur3<sizeof(_Value)>::hash((const uint8_t*) &(*first), (const uint8_t*) &(*last), seed);
    }
    
    _Value operator()(void* ptr, size_t size, _Value seed=0) const
    {
      return __dynamic_hashmurmur3<sizeof(_Value)>::hash((const uint8_t*) ptr, ((const uint8_t*) ptr) + size, seed);
    }

  };
  
};

#endif
