// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__BITPACK_DETAIL__HPP__
#define __UTILS__BITPACK_DETAIL__HPP__ 1

#include <stdint.h>

#include <algorithm>
#include <stdexcept>

namespace utils
{
  namespace bitpack
  {
    template <typename Tp, size_t ByteSize, size_t BitSize, size_t I>
    struct __struct_bitpack_loop
    {
      static const size_t unit_size = ByteSize * 8;
      static const size_t bits_total = ByteSize * BitSize * 8;
      static const size_t loop_total = ByteSize * 8;
      static const Tp     mask = (Tp(1) << BitSize) - 1;
    
      static const size_t i = loop_total - I;
      static const size_t shift_size_total = bits_total - BitSize * (i + 1);
    
      static const size_t code_pos = BitSize - 1 - (shift_size_total / unit_size);
      static const size_t shift_amount = shift_size_total & (unit_size - 1);
      static const Tp     shift_amount_mask = (Tp(1) << shift_amount) - 1;
      static const bool   filling = (unit_size - shift_amount < BitSize);
      static const bool   clearing = (unit_size - shift_amount <= BitSize);
    
      template <typename _Tp, size_t _ByteSize>
      struct __unsigned_value {};
    
      template <typename _Tp>
      struct __unsigned_value<_Tp,1> 
      {
	typedef uint8_t value_type;
      };
    
      template <typename _Tp>
      struct __unsigned_value<_Tp,2> 
      {
	typedef uint16_t value_type;
      };
    
      template <typename _Tp>
      struct __unsigned_value<_Tp,4> 
      {
	typedef uint32_t value_type;
      };
    
      template <typename _Tp>
      struct __unsigned_value<_Tp,8>
      {
	typedef uint64_t value_type;
      };
    
    
      template <typename _Tp, bool Fill>
      struct __filler {};
    
      template <typename _Tp>
      struct __filler<_Tp,true>
      {
	static inline
	void pack(const _Tp* source, _Tp* destination)
	{
	  typedef typename __unsigned_value<_Tp, ByteSize>::value_type uvalue_type;
	
	  destination[code_pos - 1] |= (uvalue_type(source[i] & mask) >> (unit_size - shift_amount)) & shift_amount_mask;
	}
      
	static inline
	void unpack(const _Tp* source, _Tp* destination)
	{
	  destination[i] = (source[code_pos - 1] << (unit_size - shift_amount)) & mask;
	}
      };
    
      template <typename _Tp>
      struct __filler<_Tp,false>
      {
	static inline
	void pack(const _Tp* source, _Tp* destination) {  }

	static inline
	void unpack(const _Tp* source, _Tp* destination) { }
      };

      template <typename _Tp, bool Equal>
      struct __assign {};
    
      template <typename _Tp>
      struct __assign<_Tp, true>
      {
	static inline
	void assign(_Tp& dest, const _Tp& value) { dest = value; }
      };
      template <typename _Tp>
      struct __assign<_Tp, false>
      {
	static inline
	void assign(_Tp& dest, const _Tp& value) { dest |= value; }
      };
    
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	__filler<Tp,filling>::pack(source, destination);
      
	__assign<Tp,clearing>::assign(destination[code_pos], (source[i] & mask) << shift_amount);
      
	__struct_bitpack_loop<Tp,ByteSize,BitSize,I-1>::pack(source, destination);
      }
    
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	typedef typename __unsigned_value<Tp, ByteSize>::value_type uvalue_type;
      
	__filler<Tp,filling>::unpack(source, destination);
      
	__assign<Tp,(! filling)>::assign(destination[i], (uvalue_type(source[code_pos]) >> shift_amount) & mask);
      
	__struct_bitpack_loop<Tp,ByteSize,BitSize,I-1>::unpack(source, destination);
      }
    };
  
    template <typename Tp, size_t ByteSize, size_t BitSize>
    struct __struct_bitpack_loop<Tp,ByteSize,BitSize,0>
    {
      static inline
      void pack(const Tp* source, Tp* destination) {}
      
      static inline
      void unpack(const Tp* source, Tp* destination) {}
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,1,1,8>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	destination[0]  = (source[0] & 0x01) << 7;
	destination[0] |= (source[1] & 0x01) << 6;
	destination[0] |= (source[2] & 0x01) << 5;
	destination[0] |= (source[3] & 0x01) << 4;
	destination[0] |= (source[4] & 0x01) << 3;
	destination[0] |= (source[5] & 0x01) << 2;
	destination[0] |= (source[6] & 0x01) << 1;
	destination[0] |= (source[7] & 0x01);
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	destination[0] = (source[0] >> 7) & 0x01;
	destination[1] = (source[0] >> 6) & 0x01;
	destination[2] = (source[0] >> 5) & 0x01;
	destination[3] = (source[0] >> 4) & 0x01;
	destination[4] = (source[0] >> 3) & 0x01;
	destination[5] = (source[0] >> 2) & 0x01;
	destination[6] = (source[0] >> 1) & 0x01;
	destination[7] = source[0] & 0x01;
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,1,2,8>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 2; ++ i) {
	  destination[i]  = (source[i * 2 + 0] & 0x03) << 6;
	  destination[i] |= (source[i * 2 + 1] & 0x03) << 4;
	  destination[i] |= (source[i * 2 + 2] & 0x03) << 2;
	  destination[i] |= (source[i * 2 + 3] & 0x03);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 2; ++ i) {
	  destination[i * 2 + 0] = (source[i] >> 6) & 0x03;
	  destination[i * 2 + 1] = (source[i] >> 4) & 0x03;
	  destination[i * 2 + 2] = (source[i] >> 2) & 0x03;
	  destination[i * 2 + 3] = source[i] & 0x03;
	}
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,1,4,8>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 4; ++ i) {
	  destination[i]  = (source[i * 2] & 0x0f) << 4;
	  destination[i] |= (source[i * 2 + 1] & 0x0f);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 4; ++ i) {
	  destination[i * 2]     = (source[i] >> 4) & 0x0f;
	  destination[i * 2 + 1] = source[i] & 0x0f;
	}
      }
    };
    
    template <typename Tp>
    struct __struct_bitpack_loop<Tp,2,1,16>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	destination[0]  = (source[0] & 0x01) << 15;
	destination[0] |= (source[1] & 0x01) << 14;
	destination[0] |= (source[2] & 0x01) << 13;
	destination[0] |= (source[3] & 0x01) << 12;
	destination[0] |= (source[4] & 0x01) << 11;
	destination[0] |= (source[5] & 0x01) << 10;
	destination[0] |= (source[6] & 0x01) << 9;
	destination[0] |= (source[7] & 0x01) << 8;
	destination[0] |= (source[8] & 0x01) << 7;
	destination[0] |= (source[9] & 0x01) << 6;
	destination[0] |= (source[10] & 0x01) << 5;
	destination[0] |= (source[11] & 0x01) << 4;
	destination[0] |= (source[12] & 0x01) << 3;
	destination[0] |= (source[13] & 0x01) << 2;
	destination[0] |= (source[14] & 0x01) << 1;
	destination[0] |= (source[15] & 0x01);
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	destination[0] = (source[0] >> 15) & 0x01;
	destination[1] = (source[0] >> 14) & 0x01;
	destination[2] = (source[0] >> 13) & 0x01;
	destination[3] = (source[0] >> 12) & 0x01;
	destination[4] = (source[0] >> 11) & 0x01;
	destination[5] = (source[0] >> 10) & 0x01;
	destination[6] = (source[0] >> 9) & 0x01;
	destination[7] = (source[0] >> 8) & 0x01;
	destination[6] = (source[0] >> 7) & 0x01;
	destination[6] = (source[0] >> 6) & 0x01;
	destination[6] = (source[0] >> 5) & 0x01;
	destination[6] = (source[0] >> 4) & 0x01;
	destination[6] = (source[0] >> 3) & 0x01;
	destination[6] = (source[0] >> 2) & 0x01;
	destination[6] = (source[0] >> 1) & 0x01;
	destination[6] = (source[0]) & 0x01;
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,2,2,16>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 2; ++ i) {
	  destination[i]  = (source[i * 8 + 0] & 0x03) << 14;
	  destination[i] |= (source[i * 8 + 1] & 0x03) << 12;
	  destination[i] |= (source[i * 8 + 2] & 0x03) << 10;
	  destination[i] |= (source[i * 8 + 3] & 0x03) << 8;
	  destination[i] |= (source[i * 8 + 4] & 0x03) << 6;
	  destination[i] |= (source[i * 8 + 5] & 0x03) << 4;
	  destination[i] |= (source[i * 8 + 6] & 0x03) << 2;
	  destination[i] |= (source[i * 8 + 7] & 0x03);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 2; ++ i) {
	  destination[i * 8 + 0] = (source[i] >> 14) & 0x03;
	  destination[i * 8 + 1] = (source[i] >> 12)  & 0x03;
	  destination[i * 8 + 2] = (source[i] >> 10)  & 0x03;
	  destination[i * 8 + 3] = (source[i] >> 8) & 0x03;
	  destination[i * 8 + 4] = (source[i] >> 6) & 0x03;
	  destination[i * 8 + 5] = (source[i] >> 4) & 0x03;
	  destination[i * 8 + 6] = (source[i] >> 2) & 0x03;
	  destination[i * 8 + 7] = source[i] & 0x03;
	}
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,2,4,16>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 4; ++ i) {
	  destination[i]  = (source[i * 4 + 0] & 0x0f) << 12;
	  destination[i] |= (source[i * 4 + 1] & 0x0f) << 8;
	  destination[i] |= (source[i * 4 + 2] & 0x0f) << 4;
	  destination[i] |= (source[i * 4 + 3] & 0x0f);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 4; ++ i) {
	  destination[i * 4 + 0] = (source[i] >> 12) & 0x0f;
	  destination[i * 4 + 1] = (source[i] >> 8)  & 0x0f;
	  destination[i * 4 + 2] = (source[i] >> 4)  & 0x0f;
	  destination[i * 4 + 3] = source[i] & 0x0f;
	}
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,2,8,16>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 8; ++ i) {
	  destination[i]  = (source[i * 2] & 0xff) << 8;
	  destination[i] |= (source[i * 2 + 1] & 0xff);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 8; ++ i) {
	  destination[i * 2]     = (source[i] >> 8) & 0xff;
	  destination[i * 2 + 1] = source[i] & 0xff;
	}
      }
    };
    
    template <typename Tp>
    struct __struct_bitpack_loop<Tp,4,1,32>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	destination[0]  = (source[0] & 0x01) << 31;
	destination[0] |= (source[1] & 0x01) << 30;
	destination[0] |= (source[2] & 0x01) << 29;
	destination[0] |= (source[3] & 0x01) << 28;
	destination[0] |= (source[4] & 0x01) << 27;
	destination[0] |= (source[5] & 0x01) << 26;
	destination[0] |= (source[6] & 0x01) << 25;
	destination[0] |= (source[7] & 0x01) << 24;
	destination[0] |= (source[8] & 0x01) << 23;
	destination[0] |= (source[9] & 0x01) << 22;
	destination[0] |= (source[10] & 0x01) << 21;
	destination[0] |= (source[11] & 0x01) << 20;
	destination[0] |= (source[12] & 0x01) << 19;
	destination[0] |= (source[13] & 0x01) << 18;
	destination[0] |= (source[14] & 0x01) << 17;
	destination[0] |= (source[15] & 0x01) << 16;
	destination[0] |= (source[16] & 0x01) << 15;
	destination[0] |= (source[17] & 0x01) << 14;
	destination[0] |= (source[18] & 0x01) << 13;
	destination[0] |= (source[19] & 0x01) << 12;
	destination[0] |= (source[20] & 0x01) << 11;
	destination[0] |= (source[21] & 0x01) << 10;
	destination[0] |= (source[22] & 0x01) << 9;
	destination[0] |= (source[23] & 0x01) << 8;
	destination[0] |= (source[24] & 0x01) << 7;
	destination[0] |= (source[25] & 0x01) << 6;
	destination[0] |= (source[26] & 0x01) << 5;
	destination[0] |= (source[27] & 0x01) << 4;
	destination[0] |= (source[28] & 0x01) << 3;
	destination[0] |= (source[29] & 0x01) << 2;
	destination[0] |= (source[30] & 0x01) << 1;
	destination[0] |= (source[31] & 0x01);
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	destination[0] = (source[0] >> 31) & 0x01;
	destination[1] = (source[0] >> 30) & 0x01;
	destination[2] = (source[0] >> 29) & 0x01;
	destination[3] = (source[0] >> 28) & 0x01;
	destination[4] = (source[0] >> 27) & 0x01;
	destination[5] = (source[0] >> 26) & 0x01;
	destination[6] = (source[0] >> 25) & 0x01;
	destination[7] = (source[0] >> 24) & 0x01;
	destination[8] = (source[0] >> 23) & 0x01;
	destination[9] = (source[0] >> 22) & 0x01;
	destination[10] = (source[0] >> 21) & 0x01;
	destination[11] = (source[0] >> 20) & 0x01;
	destination[12] = (source[0] >> 19) & 0x01;
	destination[13] = (source[0] >> 18) & 0x01;
	destination[14] = (source[0] >> 17) & 0x01;
	destination[15] = (source[0] >> 16) & 0x01;
	destination[16] = (source[0] >> 15) & 0x01;
	destination[17] = (source[0] >> 14) & 0x01;
	destination[18] = (source[0] >> 13) & 0x01;
	destination[19] = (source[0] >> 12) & 0x01;
	destination[20] = (source[0] >> 11) & 0x01;
	destination[21] = (source[0] >> 10) & 0x01;
	destination[22] = (source[0] >> 9) & 0x01;
	destination[23] = (source[0] >> 8) & 0x01;
	destination[24] = (source[0] >> 7) & 0x01;
	destination[25] = (source[0] >> 6) & 0x01;
	destination[26] = (source[0] >> 5) & 0x01;
	destination[27] = (source[0] >> 4) & 0x01;
	destination[28] = (source[0] >> 3) & 0x01;
	destination[29] = (source[0] >> 2) & 0x01;
	destination[30] = (source[0] >> 1) & 0x01;
	destination[31] = (source[0]) & 0x01;
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,4,2,32>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 2; ++ i) {
	  destination[i]  = (source[i * 16 + 0] & 0x03) << 30;
	  destination[i] |= (source[i * 16 + 1] & 0x03) << 28;
	  destination[i] |= (source[i * 16 + 2] & 0x03) << 26;
	  destination[i] |= (source[i * 16 + 3] & 0x03) << 24;
	  destination[i] |= (source[i * 16 + 4] & 0x03) << 22;
	  destination[i] |= (source[i * 16 + 5] & 0x03) << 20;
	  destination[i] |= (source[i * 16 + 6] & 0x03) << 18;
	  destination[i] |= (source[i * 16 + 7] & 0x03) << 16;
	  destination[i] |= (source[i * 16 + 8] & 0x03) << 14;
	  destination[i] |= (source[i * 16 + 9] & 0x03) << 12;
	  destination[i] |= (source[i * 16 + 10] & 0x03) << 10;
	  destination[i] |= (source[i * 16 + 11] & 0x03) << 8;
	  destination[i] |= (source[i * 16 + 12] & 0x03) << 6;
	  destination[i] |= (source[i * 16 + 13] & 0x03) << 4;
	  destination[i] |= (source[i * 16 + 14] & 0x03) << 2;
	  destination[i] |= (source[i * 16 + 15] & 0x03);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 2; ++ i) {
	  destination[i * 16 + 0] = (source[i] >> 30) & 0x03;
	  destination[i * 16 + 1] = (source[i] >> 28) & 0x03;
	  destination[i * 16 + 2] = (source[i] >> 26) & 0x03;
	  destination[i * 16 + 3] = (source[i] >> 24) & 0x03;
	  destination[i * 16 + 4] = (source[i] >> 22) & 0x03;
	  destination[i * 16 + 5] = (source[i] >> 20) & 0x03;
	  destination[i * 16 + 6] = (source[i] >> 18) & 0x03;
	  destination[i * 16 + 7] = (source[i] >> 16) & 0x03;
	  destination[i * 16 + 8] = (source[i] >> 14) & 0x03;
	  destination[i * 16 + 9] = (source[i] >> 12) & 0x03;
	  destination[i * 16 + 10] = (source[i] >> 10) & 0x03;
	  destination[i * 16 + 11] = (source[i] >> 8) & 0x03;
	  destination[i * 16 + 12] = (source[i] >> 6) & 0x03;
	  destination[i * 16 + 13] = (source[i] >> 4) & 0x03;
	  destination[i * 16 + 14] = (source[i] >> 2) & 0x03;
	  destination[i * 16 + 15] = (source[i]) & 0x03;
	}
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,4,4,32>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 4; ++ i) {
	  destination[i]  = (source[i * 8 + 0] & 0x0f) << 28;
	  destination[i] |= (source[i * 8 + 1] & 0x0f) << 24;
	  destination[i] |= (source[i * 8 + 2] & 0x0f) << 20;
	  destination[i] |= (source[i * 8 + 3] & 0x0f) << 16;
	  destination[i] |= (source[i * 8 + 4] & 0x0f) << 12;
	  destination[i] |= (source[i * 8 + 5] & 0x0f) << 8;
	  destination[i] |= (source[i * 8 + 6] & 0x0f) << 4;
	  destination[i] |= (source[i * 8 + 7] & 0x0f);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 4; ++ i) {
	  destination[i * 8 + 0] = (source[i] >> 28) & 0x0f;
	  destination[i * 8 + 1] = (source[i] >> 24) & 0x0f;
	  destination[i * 8 + 2] = (source[i] >> 20) & 0x0f;
	  destination[i * 8 + 3] = (source[i] >> 16) & 0x0f;
	  destination[i * 8 + 4] = (source[i] >> 12) & 0x0f;
	  destination[i * 8 + 5] = (source[i] >> 8) & 0x0f;
	  destination[i * 8 + 6] = (source[i] >> 4) & 0x0f;
	  destination[i * 8 + 7] = (source[i]) & 0x0f;
	}
      }
    };
    
    template <typename Tp>
    struct __struct_bitpack_loop<Tp,4,8,32>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 8; ++ i) {
	  destination[i]  = (source[i * 4 + 0] & 0xff) << 24;
	  destination[i] |= (source[i * 4 + 1] & 0xff) << 16;
	  destination[i] |= (source[i * 4 + 2] & 0xff) << 8;
	  destination[i] |= (source[i * 4 + 3] & 0xff);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 8; ++ i) {
	  destination[i * 4 + 0] = (source[i] >> 24) & 0xff;
	  destination[i * 4 + 1] = (source[i] >> 16) & 0xff;
	  destination[i * 4 + 2] = (source[i] >> 8) & 0xff;
	  destination[i * 4 + 3] = source[i] & 0xff;
	}
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,4,16,32>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 16; ++ i) {
	  destination[i]  = (source[i * 2] & 0xffff) << 16;
	  destination[i] |= (source[i * 2 + 1] & 0xffff);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 16; ++ i) {
	  destination[i * 2]     = (source[i] >> 16) & 0xffff;
	  destination[i * 2 + 1] = source[i] & 0xffff;
	}
      }
    };
    
    template <typename Tp>
    struct __struct_bitpack_loop<Tp,8,1,64>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	destination[0]  = (source[0] & 0x01) << 63;
	destination[0] |= (source[1] & 0x01) << 62;
	destination[0] |= (source[2] & 0x01) << 61;
	destination[0] |= (source[3] & 0x01) << 60;
	destination[0] |= (source[4] & 0x01) << 59;
	destination[0] |= (source[5] & 0x01) << 58;
	destination[0] |= (source[6] & 0x01) << 57;
	destination[0] |= (source[7] & 0x01) << 56;
	destination[0] |= (source[8] & 0x01) << 55;
	destination[0] |= (source[9] & 0x01) << 54;
	destination[0] |= (source[10] & 0x01) << 53;
	destination[0] |= (source[11] & 0x01) << 52;
	destination[0] |= (source[12] & 0x01) << 51;
	destination[0] |= (source[13] & 0x01) << 50;
	destination[0] |= (source[14] & 0x01) << 49;
	destination[0] |= (source[15] & 0x01) << 48;
	destination[0] |= (source[16] & 0x01) << 47;
	destination[0] |= (source[17] & 0x01) << 46;
	destination[0] |= (source[18] & 0x01) << 45;
	destination[0] |= (source[19] & 0x01) << 44;
	destination[0] |= (source[20] & 0x01) << 43;
	destination[0] |= (source[21] & 0x01) << 42;
	destination[0] |= (source[22] & 0x01) << 41;
	destination[0] |= (source[23] & 0x01) << 40;
	destination[0] |= (source[24] & 0x01) << 39;
	destination[0] |= (source[25] & 0x01) << 38;
	destination[0] |= (source[26] & 0x01) << 37;
	destination[0] |= (source[27] & 0x01) << 36;
	destination[0] |= (source[28] & 0x01) << 35;
	destination[0] |= (source[29] & 0x01) << 34;
	destination[0] |= (source[30] & 0x01) << 33;
	destination[0] |= (source[31] & 0x01) << 32;
	destination[0] |= (source[32] & 0x01) << 31;
	destination[0] |= (source[33] & 0x01) << 30;
	destination[0] |= (source[34] & 0x01) << 29;
	destination[0] |= (source[35] & 0x01) << 28;
	destination[0] |= (source[36] & 0x01) << 27;
	destination[0] |= (source[37] & 0x01) << 26;
	destination[0] |= (source[38] & 0x01) << 25;
	destination[0] |= (source[39] & 0x01) << 24;
	destination[0] |= (source[40] & 0x01) << 23;
	destination[0] |= (source[41] & 0x01) << 22;
	destination[0] |= (source[42] & 0x01) << 21;
	destination[0] |= (source[43] & 0x01) << 20;
	destination[0] |= (source[44] & 0x01) << 19;
	destination[0] |= (source[45] & 0x01) << 18;
	destination[0] |= (source[46] & 0x01) << 17;
	destination[0] |= (source[47] & 0x01) << 16;
	destination[0] |= (source[48] & 0x01) << 15;
	destination[0] |= (source[49] & 0x01) << 14;
	destination[0] |= (source[50] & 0x01) << 13;
	destination[0] |= (source[51] & 0x01) << 12;
	destination[0] |= (source[52] & 0x01) << 11;
	destination[0] |= (source[53] & 0x01) << 10;
	destination[0] |= (source[54] & 0x01) << 9;
	destination[0] |= (source[55] & 0x01) << 8;
	destination[0] |= (source[56] & 0x01) << 7;
	destination[0] |= (source[57] & 0x01) << 6;
	destination[0] |= (source[58] & 0x01) << 5;
	destination[0] |= (source[59] & 0x01) << 4;
	destination[0] |= (source[60] & 0x01) << 3;
	destination[0] |= (source[61] & 0x01) << 2;
	destination[0] |= (source[62] & 0x01) << 1;
	destination[0] |= (source[63] & 0x01);
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	destination[0] = (source[0] >> 63) & 0x01;
	destination[1] = (source[0] >> 62) & 0x01;
	destination[2] = (source[0] >> 61) & 0x01;
	destination[3] = (source[0] >> 60) & 0x01;
	destination[4] = (source[0] >> 59) & 0x01;
	destination[5] = (source[0] >> 58) & 0x01;
	destination[6] = (source[0] >> 57) & 0x01;
	destination[7] = (source[0] >> 56) & 0x01;
	destination[8] = (source[0] >> 55) & 0x01;
	destination[9] = (source[0] >> 54) & 0x01;
	destination[10] = (source[0] >> 53) & 0x01;
	destination[11] = (source[0] >> 52) & 0x01;
	destination[12] = (source[0] >> 51) & 0x01;
	destination[13] = (source[0] >> 50) & 0x01;
	destination[14] = (source[0] >> 49) & 0x01;
	destination[15] = (source[0] >> 48) & 0x01;
	destination[16] = (source[0] >> 47) & 0x01;
	destination[17] = (source[0] >> 46) & 0x01;
	destination[18] = (source[0] >> 45) & 0x01;
	destination[19] = (source[0] >> 44) & 0x01;
	destination[20] = (source[0] >> 43) & 0x01;
	destination[21] = (source[0] >> 42) & 0x01;
	destination[22] = (source[0] >> 41) & 0x01;
	destination[23] = (source[0] >> 40) & 0x01;
	destination[24] = (source[0] >> 39) & 0x01;
	destination[25] = (source[0] >> 38) & 0x01;
	destination[26] = (source[0] >> 37) & 0x01;
	destination[27] = (source[0] >> 36) & 0x01;
	destination[28] = (source[0] >> 35) & 0x01;
	destination[29] = (source[0] >> 34) & 0x01;
	destination[30] = (source[0] >> 33) & 0x01;
	destination[31] = (source[0] >> 32) & 0x01;
	destination[32] = (source[0] >> 31) & 0x01;
	destination[33] = (source[0] >> 30) & 0x01;
	destination[34] = (source[0] >> 29) & 0x01;
	destination[35] = (source[0] >> 28) & 0x01;
	destination[36] = (source[0] >> 27) & 0x01;
	destination[37] = (source[0] >> 26) & 0x01;
	destination[38] = (source[0] >> 25) & 0x01;
	destination[39] = (source[0] >> 24) & 0x01;
	destination[40] = (source[0] >> 23) & 0x01;
	destination[41] = (source[0] >> 22) & 0x01;
	destination[42] = (source[0] >> 21) & 0x01;
	destination[43] = (source[0] >> 20) & 0x01;
	destination[44] = (source[0] >> 19) & 0x01;
	destination[45] = (source[0] >> 18) & 0x01;
	destination[46] = (source[0] >> 17) & 0x01;
	destination[47] = (source[0] >> 16) & 0x01;
	destination[48] = (source[0] >> 15) & 0x01;
	destination[49] = (source[0] >> 14) & 0x01;
	destination[50] = (source[0] >> 13) & 0x01;
	destination[51] = (source[0] >> 12) & 0x01;
	destination[52] = (source[0] >> 11) & 0x01;
	destination[53] = (source[0] >> 10) & 0x01;
	destination[54] = (source[0] >> 9) & 0x01;
	destination[55] = (source[0] >> 8) & 0x01;
	destination[56] = (source[0] >> 7) & 0x01;
	destination[57] = (source[0] >> 6) & 0x01;
	destination[58] = (source[0] >> 5) & 0x01;
	destination[59] = (source[0] >> 4) & 0x01;
	destination[60] = (source[0] >> 3) & 0x01;
	destination[61] = (source[0] >> 2) & 0x01;
	destination[62] = (source[0] >> 1) & 0x01;
	destination[63] = (source[0]) & 0x01;
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,8,2,64>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 2; ++ i) {
	  destination[i]  = (source[i * 32 + 0] & 0x03) << 62;
	  destination[i] |= (source[i * 32 + 1] & 0x03) << 60;
	  destination[i] |= (source[i * 32 + 2] & 0x03) << 58;
	  destination[i] |= (source[i * 32 + 3] & 0x03) << 56;
	  destination[i] |= (source[i * 32 + 4] & 0x03) << 54;
	  destination[i] |= (source[i * 32 + 5] & 0x03) << 52;
	  destination[i] |= (source[i * 32 + 6] & 0x03) << 50;
	  destination[i] |= (source[i * 32 + 7] & 0x03) << 48;
	  destination[i] |= (source[i * 32 + 8] & 0x03) << 46;
	  destination[i] |= (source[i * 32 + 9] & 0x03) << 44;
	  destination[i] |= (source[i * 32 + 10] & 0x03) << 42;
	  destination[i] |= (source[i * 32 + 11] & 0x03) << 40;
	  destination[i] |= (source[i * 32 + 12] & 0x03) << 38;
	  destination[i] |= (source[i * 32 + 13] & 0x03) << 36;
	  destination[i] |= (source[i * 32 + 14] & 0x03) << 34;
	  destination[i] |= (source[i * 32 + 15] & 0x03) << 32;
	  destination[i] |= (source[i * 32 + 16] & 0x03) << 30;
	  destination[i] |= (source[i * 32 + 17] & 0x03) << 28;
	  destination[i] |= (source[i * 32 + 18] & 0x03) << 26;
	  destination[i] |= (source[i * 32 + 19] & 0x03) << 24;
	  destination[i] |= (source[i * 32 + 20] & 0x03) << 22;
	  destination[i] |= (source[i * 32 + 21] & 0x03) << 20;
	  destination[i] |= (source[i * 32 + 22] & 0x03) << 18;
	  destination[i] |= (source[i * 32 + 23] & 0x03) << 16;
	  destination[i] |= (source[i * 32 + 24] & 0x03) << 14;
	  destination[i] |= (source[i * 32 + 25] & 0x03) << 12;
	  destination[i] |= (source[i * 32 + 26] & 0x03) << 10;
	  destination[i] |= (source[i * 32 + 27] & 0x03) << 8;
	  destination[i] |= (source[i * 32 + 28] & 0x03) << 6;
	  destination[i] |= (source[i * 32 + 29] & 0x03) << 4;
	  destination[i] |= (source[i * 32 + 30] & 0x03) << 2;
	  destination[i] |= (source[i * 32 + 31] & 0x03);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 2; ++ i) {
	  destination[i * 32 + 0] = (source[i] >> 62) & 0x03;
	  destination[i * 32 + 1] = (source[i] >> 60) & 0x03;
	  destination[i * 32 + 2] = (source[i] >> 58) & 0x03;
	  destination[i * 32 + 3] = (source[i] >> 56) & 0x03;
	  destination[i * 32 + 4] = (source[i] >> 54) & 0x03;
	  destination[i * 32 + 5] = (source[i] >> 52) & 0x03;
	  destination[i * 32 + 6] = (source[i] >> 50) & 0x03;
	  destination[i * 32 + 7] = (source[i] >> 48) & 0x03;
	  destination[i * 32 + 8] = (source[i] >> 46) & 0x03;
	  destination[i * 32 + 9] = (source[i] >> 44) & 0x03;
	  destination[i * 32 + 10] = (source[i] >> 42) & 0x03;
	  destination[i * 32 + 11] = (source[i] >> 40) & 0x03;
	  destination[i * 32 + 12] = (source[i] >> 38) & 0x03;
	  destination[i * 32 + 13] = (source[i] >> 36) & 0x03;
	  destination[i * 32 + 14] = (source[i] >> 34) & 0x03;
	  destination[i * 32 + 15] = (source[i] >> 32) & 0x03;
	  destination[i * 32 + 16] = (source[i] >> 30) & 0x03;
	  destination[i * 32 + 17] = (source[i] >> 28) & 0x03;
	  destination[i * 32 + 18] = (source[i] >> 26) & 0x03;
	  destination[i * 32 + 19] = (source[i] >> 24) & 0x03;
	  destination[i * 32 + 20] = (source[i] >> 22) & 0x03;
	  destination[i * 32 + 21] = (source[i] >> 20) & 0x03;
	  destination[i * 32 + 22] = (source[i] >> 18) & 0x03;
	  destination[i * 32 + 23] = (source[i] >> 16) & 0x03;
	  destination[i * 32 + 24] = (source[i] >> 14) & 0x03;
	  destination[i * 32 + 25] = (source[i] >> 12) & 0x03;
	  destination[i * 32 + 26] = (source[i] >> 10) & 0x03;
	  destination[i * 32 + 27] = (source[i] >> 8) & 0x03;
	  destination[i * 32 + 28] = (source[i] >> 6) & 0x03;
	  destination[i * 32 + 29] = (source[i] >> 4) & 0x03;
	  destination[i * 32 + 30] = (source[i] >> 2) & 0x03;
	  destination[i * 32 + 31] = (source[i]) & 0x03;
	}
      }
    };
    
    template <typename Tp>
    struct __struct_bitpack_loop<Tp,8,4,64>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 4; ++ i) {
	  destination[i]  = (source[i * 16 + 0] & 0x0f) << 60;
	  destination[i] |= (source[i * 16 + 1] & 0x0f) << 56;
	  destination[i] |= (source[i * 16 + 2] & 0x0f) << 52;
	  destination[i] |= (source[i * 16 + 3] & 0x0f) << 48;
	  destination[i] |= (source[i * 16 + 4] & 0x0f) << 44;
	  destination[i] |= (source[i * 16 + 5] & 0x0f) << 40;
	  destination[i] |= (source[i * 16 + 6] & 0x0f) << 36;
	  destination[i] |= (source[i * 16 + 7] & 0x0f) << 32;
	  destination[i] |= (source[i * 16 + 8] & 0x0f) << 28;
	  destination[i] |= (source[i * 16 + 9] & 0x0f) << 24;
	  destination[i] |= (source[i * 16 + 10] & 0x0f) << 20;
	  destination[i] |= (source[i * 16 + 11] & 0x0f) << 16;
	  destination[i] |= (source[i * 16 + 12] & 0x0f) << 12;
	  destination[i] |= (source[i * 16 + 13] & 0x0f) << 8;
	  destination[i] |= (source[i * 16 + 14] & 0x0f) << 4;
	  destination[i] |= (source[i * 16 + 15] & 0x0f);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 4; ++ i) {
	  destination[i * 16 + 0] = (source[i] >> 60) & 0x0f;
	  destination[i * 16 + 1] = (source[i] >> 56) & 0x0f;
	  destination[i * 16 + 2] = (source[i] >> 52) & 0x0f;
	  destination[i * 16 + 3] = (source[i] >> 48) & 0x0f;
	  destination[i * 16 + 4] = (source[i] >> 44) & 0x0f;
	  destination[i * 16 + 5] = (source[i] >> 40) & 0x0f;
	  destination[i * 16 + 6] = (source[i] >> 36) & 0x0f;
	  destination[i * 16 + 7] = (source[i] >> 32) & 0x0f;
	  destination[i * 16 + 8] = (source[i] >> 28) & 0x0f;
	  destination[i * 16 + 9] = (source[i] >> 24) & 0x0f;
	  destination[i * 16 + 10] = (source[i] >> 20) & 0x0f;
	  destination[i * 16 + 11] = (source[i] >> 16) & 0x0f;
	  destination[i * 16 + 12] = (source[i] >> 12) & 0x0f;
	  destination[i * 16 + 13] = (source[i] >> 8) & 0x0f;
	  destination[i * 16 + 14] = (source[i] >> 4) & 0x0f;
	  destination[i * 16 + 15] = (source[i]) & 0x0f;
	}
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,8,8,64>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 8; ++ i) {
	  destination[i]  = (source[i * 8 + 0] & 0xff) << 56;
	  destination[i] |= (source[i * 8 + 1] & 0xff) << 48;
	  destination[i] |= (source[i * 8 + 2] & 0xff) << 40;
	  destination[i] |= (source[i * 8 + 3] & 0xff) << 32;
	  destination[i] |= (source[i * 8 + 4] & 0xff) << 24;
	  destination[i] |= (source[i * 8 + 5] & 0xff) << 16;
	  destination[i] |= (source[i * 8 + 6] & 0xff) << 8;
	  destination[i] |= (source[i * 8 + 7] & 0xff);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 8; ++ i) {
	  destination[i * 8 + 0] = (source[i] >> 56) & 0xff;
	  destination[i * 8 + 1] = (source[i] >> 48) & 0xff;
	  destination[i * 8 + 2] = (source[i] >> 40) & 0xff;
	  destination[i * 8 + 3] = (source[i] >> 32) & 0xff;
	  destination[i * 8 + 4] = (source[i] >> 24) & 0xff;
	  destination[i * 8 + 5] = (source[i] >> 16) & 0xff;
	  destination[i * 8 + 6] = (source[i] >> 8) & 0xff;
	  destination[i * 8 + 7] = (source[i]) & 0xff;
	}
      }
    };

    template <typename Tp>
    struct __struct_bitpack_loop<Tp,8,16,64>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 16; ++ i) {
	  destination[i]  = (source[i * 4 + 0] & 0xffff) << 48;
	  destination[i] |= (source[i * 4 + 1] & 0xffff) << 32;
	  destination[i] |= (source[i * 4 + 2] & 0xffff) << 16;
	  destination[i] |= (source[i * 4 + 3] & 0xffff);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 16; ++ i) {
	  destination[i * 4 + 0] = (source[i] >> 48) & 0xffff;
	  destination[i * 4 + 1] = (source[i] >> 32) & 0xffff;
	  destination[i * 4 + 2] = (source[i] >> 16) & 0xffff;
	  destination[i * 4 + 3] = (source[i]) & 0xffff;
	}
      }
    };
    
    template <typename Tp>
    struct __struct_bitpack_loop<Tp,8,32,64>
    {
      static inline
      void pack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 32; ++ i) {
	  destination[i]  = (source[i * 2] & 0xffffffffull) << 32;
	  destination[i] |= (source[i * 2 + 1] & 0xffffffffull);
	}
      }
      
      static inline
      void unpack(const Tp* source, Tp* destination)
      {
	for (int i = 0; i < 32; ++ i) {
	  destination[i * 2]     = (source[i] >> 32) & 0xffffffffull;
	  destination[i * 2 + 1] = source[i] & 0xffffffffull;
	}
      }
    };
  };
};

#endif
