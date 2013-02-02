// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//
// this is a reimplementation of xxHash:
/*
   xxHash - Fast Hash algorithm
   Header File
   Copyright (C) 2012, Yann Collet.
   BSD 2-Clause License (http://www.opensource.org/licenses/bsd-license.php)

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
  
       * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
       * Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the
   distribution.
  
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

        You can contact the author at :
        - xxHash source repository : http://code.google.com/p/xxhash/
*/


#ifndef __UTILS__HASHXX__HPP__
#define __UTILS__HASHXX__HPP__ 1

#include <cstddef>
#include <stdint.h>

namespace utils
{
  
  struct __base_hashxx
  {
    static inline
    uint32_t XXH_rotl32(uint32_t x, uint32_t r)
    {
      return (x << r) | (x >> (32 - r));
    }

    template <typename Iterator>
    static inline
    uint32_t XXH_LE32(Iterator p)
    {
      return *((uint32_t*) p);
    }

    static const uint32_t PRIME32_1 = 2654435761U;
    static const uint32_t PRIME32_2 = 2246822519U;
    static const uint32_t PRIME32_3 = 3266489917U;
    static const uint32_t PRIME32_4 =  668265263U;
    static const uint32_t PRIME32_5 =  374761393U;
  };

  template <size_t Loop>
  struct __static_hashxx16 : __base_hashxx
  {
    template <typename Iterator>
    static inline
    uint32_t hash(Iterator& p, uint32_t seed)
    {
      uint32_t v1 = seed + PRIME32_1 + PRIME32_2;
      uint32_t v2 = seed + PRIME32_2;
      uint32_t v3 = seed + 0;
      uint32_t v4 = seed - PRIME32_1;
      
      for (size_t i = 0; i != Loop; ++ i) {
	v1 += XXH_LE32(p)      * PRIME32_2; v1 = XXH_rotl32(v1, 13); v1 *= PRIME32_1;
	v2 += XXH_LE32(p +  4) * PRIME32_2; v2 = XXH_rotl32(v2, 13); v2 *= PRIME32_1;
	v3 += XXH_LE32(p +  8) * PRIME32_2; v3 = XXH_rotl32(v3, 13); v3 *= PRIME32_1;
	v4 += XXH_LE32(p + 12) * PRIME32_2; v4 = XXH_rotl32(v4, 13); v4 *= PRIME32_1;
	p += 16;
      };
      
      return XXH_rotl32(v1, 1) + XXH_rotl32(v2, 7) + XXH_rotl32(v3, 12) + XXH_rotl32(v4, 18);
    }
  };
  
  template <>
  struct __static_hashxx16<0> : __base_hashxx
  {
    template <typename Iterator>
    static inline
    uint32_t hash(Iterator& p, uint32_t seed)
    {
      return seed + PRIME32_5;
    }
  };

  template <size_t Loop, size_t Remain>
  struct __static_hashxx4 : __base_hashxx
  {
    template <typename Iterator>
    static inline
    void hash(Iterator& p, uint32_t& h32)
    {
      h32 += XXH_LE32(p + ((Remain - Loop) << 2)) * PRIME32_3;
      h32 = XXH_rotl32(h32, 17) * PRIME32_4;
      
      __static_hashxx4<Loop - 1, Remain>::hash(p, h32);
    }
  };

  template <size_t Remain>
  struct __static_hashxx4<0, Remain>
  {
    template <typename Iterator>
    static inline
    void hash(Iterator& p, uint32_t& h32)
    {
      p += (Remain << 2);
    }
  };
  
  template <size_t Loop, size_t Remain>
  struct __static_hashxx1 : __base_hashxx
  {
    template <typename Iterator>
    static inline
    void hash(Iterator& p, uint32_t& h32)
    {
      h32 += (*(p + Remain - Loop)) * PRIME32_5;
      h32 = XXH_rotl32(h32, 11) * PRIME32_1 ;
      
      __static_hashxx1<Loop - 1, Remain>::hash(p, h32);
    }
  };

  template <size_t Remain>
  struct __static_hashxx1<0, Remain>
  {
    template <typename Iterator>
    static inline
    void hash(Iterator& p, uint32_t& h32)
    {
      p += Remain;
    }
  };
  
  template <size_t Len>
  struct __static_hashxx_impl : __base_hashxx
  {
    static inline
    uint32_t hash(const uint8_t* p, uint32_t seed)
    {
      uint32_t h32 = __static_hashxx16<(Len >> 4)>::hash(p, seed) + static_cast<uint32_t>(Len);
      
      __static_hashxx4<((Len >> 2) & 0x03), ((Len >> 2) & 0x03)>::hash(p, h32);
      
      __static_hashxx1<(Len & 0x03), (Len & 0x03)>::hash(p, h32);
      
      h32 ^= h32 >> 15;
      h32 *= PRIME32_2;
      h32 ^= h32 >> 13;
      h32 *= PRIME32_3;
      h32 ^= h32 >> 16;
      
      return h32; 
    }
  };
  
  template <size_t _Length, size_t _Size>
  struct __static_hashxx {};

  template <size_t _Length>
  struct __static_hashxx<_Length, 4>
  {
    static inline
    uint32_t hash(const uint8_t* data, uint32_t seed)
    {
      return __static_hashxx_impl<_Length>::hash(data, seed);
    }
  };
  
  template <size_t _Length>
  struct __static_hashxx<_Length, 8>
  {
    static inline
    uint64_t hash(const uint8_t* data, uint64_t seed)
    {
      return __static_hashxx_impl<_Length>::hash(data, seed);
    }
  };
  
  template <size_t _Size>
  struct __dynamic_hashxx : __base_hashxx {};

  template <>
  struct __dynamic_hashxx<4> : __base_hashxx
  {
    template <typename Iterator>
    static inline
    uint32_t hash(Iterator p, Iterator bEnd, uint32_t seed)
    {
      const size_t len = bEnd - p;
      uint32_t h32;
      
      if (len >= 16) {
	const uint8_t* const limit = bEnd - 16;
	uint32_t v1 = seed + PRIME32_1 + PRIME32_2;
	uint32_t v2 = seed + PRIME32_2;
	uint32_t v3 = seed + 0;
	uint32_t v4 = seed - PRIME32_1;
	
	do {
	  v1 += XXH_LE32(p)      * PRIME32_2; v1 = XXH_rotl32(v1, 13); v1 *= PRIME32_1;
	  v2 += XXH_LE32(p +  4) * PRIME32_2; v2 = XXH_rotl32(v2, 13); v2 *= PRIME32_1;
	  v3 += XXH_LE32(p +  8) * PRIME32_2; v3 = XXH_rotl32(v3, 13); v3 *= PRIME32_1;
	  v4 += XXH_LE32(p + 12) * PRIME32_2; v4 = XXH_rotl32(v4, 13); v4 *= PRIME32_1;
	  p += 16;
	} while (p <= limit) ;
	
	h32 = XXH_rotl32(v1, 1) + XXH_rotl32(v2, 7) + XXH_rotl32(v3, 12) + XXH_rotl32(v4, 18);
      } else
	h32  = seed + PRIME32_5;
      
      h32 += static_cast<uint32_t>(len);

      const size_t remain4 = (len >> 2) & 0x03;
      switch (remain4) {
      case 3:
	h32 += XXH_LE32(p + ((remain4 - 3) << 2)) * PRIME32_3;
	h32 = XXH_rotl32(h32, 17) * PRIME32_4;
      case 2:
	h32 += XXH_LE32(p + ((remain4 - 2) << 2)) * PRIME32_3;
	h32 = XXH_rotl32(h32, 17) * PRIME32_4;
      case 1:
	h32 += XXH_LE32(p + ((remain4 - 1) << 2)) * PRIME32_3;
	h32 = XXH_rotl32(h32, 17) * PRIME32_4;
      }

      p += (remain4 << 2);
      
      const size_t remain1 = len & 0x03;
      switch (remain1) {
      case 3:
	h32 += (*(p + remain1 - 3)) * PRIME32_5;
	h32 = XXH_rotl32(h32, 11) * PRIME32_1;
      case 2:
	h32 += (*(p + remain1 - 2)) * PRIME32_5;
	h32 = XXH_rotl32(h32, 11) * PRIME32_1;
      case 1:
	h32 += (*(p + remain1 - 1)) * PRIME32_5;
	h32 = XXH_rotl32(h32, 11) * PRIME32_1;
      }

      p += remain1;
      
      h32 ^= h32 >> 15;
      h32 *= PRIME32_2;
      h32 ^= h32 >> 13;
      h32 *= PRIME32_3;
      h32 ^= h32 >> 16;

      return h32; 
    }
  };
  
  template <>
  struct __dynamic_hashxx<8>
  {
    template <typename Iterator>
    static inline
    uint64_t hash(Iterator p, Iterator bEnd, uint64_t seed)
    {
      return __dynamic_hashxx<4>::hash(p, bEnd, seed);
    }
  };
  
  
  
  template <typename _Value>
  struct hashxx
  {
    
  public:
    template <typename Tp, size_t N>
    _Value operator()(const Tp (&x)[N], _Value seed=0) const
    { 
      return __static_hashxx<sizeof(Tp) * N, sizeof(_Value)>::hash((const uint8_t*) x, seed);
    }
    
    template <typename Tp>
    _Value operator()(const Tp& x, _Value seed=0) const
    {
      return __static_hashxx<sizeof(Tp), sizeof(_Value)>::hash((const uint8_t*) &x, seed);
    }
    
    template <typename Iterator>
    _Value operator()(Iterator first, Iterator last, _Value seed) const
    {
      return __dynamic_hashxx<sizeof(_Value)>::hash((const uint8_t*) &(*first), (const uint8_t*) &(*last), seed);
    }
    
    _Value operator()(void* ptr, size_t size, _Value seed=0) const
    {
      return __dynamic_hashxx<sizeof(_Value)>::hash((const uint8_t*) ptr, ((const uint8_t*) ptr) + size, seed);
    }

  };
};

#endif
