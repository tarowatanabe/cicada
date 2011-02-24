// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__BITHACK__HPP__
#define __UTILS__BITHACK__HPP__

#include <stdint.h>

namespace utils
{
  namespace bithack
  {
    template <typename Tp>
    inline
    Tp branch(const bool cond, const Tp& x, const Tp& y)
    {
      const Tp mask = Tp(cond) - 1;
      return ((~mask) & x) | (mask & y);
    }
    
    template <typename Tp>
    inline
    Tp abs(Tp x)
    {
      const Tp mask = x >> (sizeof(Tp) * 8 - 1);
      return (x ^ mask) - mask;
    }

    template <typename Tp>
    inline
    Tp average(Tp x, Tp y)
    {
      return (x & y) + ((x ^ y) >> 1);
    }

    template <typename Tp>
    inline
    Tp max(Tp x, Tp y)
    {
      // integral only max..
      return x - ((x - y) & -(x < y));
    }

    template <typename Tp>
    inline
    Tp min(Tp x, Tp y)
    {
      // integral only min...
      return y + ((x - y) & -(x < y));
    }

    // we will implement metaprogramming and conventional version
    
    template <size_t ByteSize>
    struct struct_is_power2
    {
      template <typename Tp>
      bool operator()(Tp x) const { return ! (x & (x - 1)); }
    };
    
    // is-power2
    template <typename Tp>
    inline 
    bool is_power2(Tp x) { return ! (x & (x - 1)); }
    
    template <uint64_t X>
    struct static_is_power2
    {
      static const bool result = ! (X & (X - 1));
      static const bool value = ! (X & (X - 1));
    };
    
    // next-largest-power2
    template <size_t ByteSize>
    struct struct_next_largest_power2 {};
    
    template <>
    struct struct_next_largest_power2 <1>
    {
      uint8_t operator()(uint8_t x) const
      {
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	return x + 1;
      }
    };
    
    template <>
    struct struct_next_largest_power2 <2>
    {
      uint16_t operator()(uint16_t x) const
      {
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	return x + 1;
      }
    };
    
    template <>
    struct struct_next_largest_power2 <4>
    {
      uint32_t operator()(uint32_t x) const
      {
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
	return x + 1;
      }
    };
    
    template <>
    struct struct_next_largest_power2 <8>
    {
      uint64_t operator()(uint64_t x) const
      {
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
	x |= (x >> 32);
	return x + 1;
      }
    };
    
    template <typename Tp>
    inline
    uint64_t next_largest_power2(Tp x)
    {
      static struct_next_largest_power2<sizeof(Tp)> __func;
      return __func(x);
    }
    
    template <uint64_t X>
    struct static_next_largest_power2
    {
    private:
      static const uint64_t X1 = uint64_t(X)  | (uint64_t(X)  >>  1);
      static const uint64_t X2 = uint64_t(X1) | (uint64_t(X1) >>  2);
      static const uint64_t X3 = uint64_t(X2) | (uint64_t(X2) >>  4);
      static const uint64_t X4 = uint64_t(X3) | (uint64_t(X3) >>  8);
      static const uint64_t X5 = uint64_t(X4) | (uint64_t(X4) >> 16);
      static const uint64_t X6 = uint64_t(X5) | (uint64_t(X5) >> 32);
      
    public:
      static const uint64_t result = uint64_t(X6 + 1);
      static const uint64_t value = uint64_t(X6 + 1);
    };

    // most significant bit...
    template <size_t ByteCount>
    struct struct_most_significant_bit {};
    
    template <>
    struct struct_most_significant_bit <1>
    {
      uint8_t operator()(uint8_t x) const
      {
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	return (x & ~(x >> 1));
      }
    };
    
    template <>
    struct struct_most_significant_bit <2>
    {
      uint16_t operator()(uint16_t x) const
      {
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	return (x & ~(x >> 1));
      }
    };
    
    template <>
    struct struct_most_significant_bit <4>
    {
      uint32_t operator()(uint32_t x) const
      {
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
	return (x & ~(x >> 1));
      }
    };
    
    template <>
    struct struct_most_significant_bit <8>
    {
      uint64_t operator()(uint64_t x) const
      {
	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
	x |= (x >> 32);
	return (x & ~(x >> 1));
      }
    };
    
    template <typename Tp>
    inline
    uint64_t most_significant_bit(Tp x)
    {
      static struct_most_significant_bit<sizeof(Tp)> __func;
      return __func(x);
    }
    

    template <uint64_t X>
    struct static_most_significant_bit
    {
    private:
      static const uint64_t X1 = uint64_t(X)  | (uint64_t(X)  >>  1);
      static const uint64_t X2 = uint64_t(X1) | (uint64_t(X1) >>  2);
      static const uint64_t X3 = uint64_t(X2) | (uint64_t(X2) >>  4);
      static const uint64_t X4 = uint64_t(X3) | (uint64_t(X3) >>  8);
      static const uint64_t X5 = uint64_t(X4) | (uint64_t(X4) >> 16);
      static const uint64_t X6 = uint64_t(X5) | (uint64_t(X5) >> 32);
    public:
      static const uint64_t result = (X6 & ~(X6 >> 1));
      static const uint64_t value = (X6 & ~(X6 >> 1));
    };
    
    // bit-count
    // for 8-bit and 16-bit, use table-lookup
    // for 32-bit and 64-bit, use pos-count via masking/shifting etc...

    struct __bit_count_mask
    {
      const uint8_t* operator()() const
      {
	static uint8_t mask_blocks[256] = {
	  0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
	  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
	  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
	  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
	  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
	  1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
	  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
	  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
	  3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
	  4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
	};
	return mask_blocks;
      }
    };
    
    template <size_t ByteSize>
    struct struct_bit_count {};
    
    template <>
    struct struct_bit_count<1>
    {
      size_t operator()(uint8_t x) const
      {
	return static_cast<size_t>(__mask()[x]);
      }
      
      __bit_count_mask __mask;
    };
    
    template <>
    struct struct_bit_count<2>
    {
      size_t operator()(uint16_t x) const
      {
	return __func((x >> 8) & 0xff) + __func(x & 0xff);
      }
      
      struct_bit_count<1> __func;
    };
    
    template <>
    struct struct_bit_count<4>
    {
      size_t operator()(uint32_t x) const
      {
	// from http://graphics.stanford.edu/~seander/bithacks.html
	
	x = x - ((x >> 1) & 0x55555555);                    // reuse input as temporary
	x = (x & 0x33333333) + ((x >> 2) & 0x33333333);     // temp
	return static_cast<size_t>((((x + (x >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24);
      }
    };
    
    template <>
    struct struct_bit_count<8>
    {
      size_t operator()(uint64_t x) const
      {
	x = x - ((x >> 1) & 0x5555555555555555LLU);
	x = (x & 0x3333333333333333LLU) + ((x >> 2) & 0x3333333333333333LLU);
	x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FLLU;
	x = x + (x >> 8);
	x = x + (x >> 16);
	x = x + (x >> 32);
	
	return static_cast<size_t>(x & 0xFF);
      }
    };
    
    template <typename Tp>
    inline
    size_t bit_count(Tp x)
    {
      static struct_bit_count<sizeof(Tp)> __func;
      return __func(x);
    }
    
    template <uint64_t X>
    struct static_bit_count
    {
    private:
      static const uint64_t X1 = X - ((X >> 1) & 0x5555555555555555LLU);
      static const uint64_t X2 = (X1 & 0x3333333333333333LLU) + ((X1 >> 2) & 0x3333333333333333LLU);
      static const uint64_t X3 = (X2 + (X2 >> 4)) & 0x0F0F0F0F0F0F0F0FLLU;
      static const uint64_t X4 = X3 + (X3 >> 8);
      static const uint64_t X5 = X4 + (X4 >> 16);
      static const uint64_t X6 = X5 + (X5 >> 32);
    public:
      static const uint64_t result = X6 & 0xFF;
      static const uint64_t value  = X6 & 0xFF;
    };
    
    // floor_log2
    template <size_t ByteSize>
    struct struct_floor_log2 {};
    
    template <>
    struct struct_floor_log2<1>
    {
      uint8_t operator()(uint8_t x) const
      {
	x |= (x >> 1);
        x |= (x >> 2);
        x |= (x >> 4);
	return (bit_count(x >> 1));
      }
    };
    
    template <>
    struct struct_floor_log2<2>
    {
      uint16_t operator()(uint16_t x) const
      {
	x |= (x >> 1);
        x |= (x >> 2);
        x |= (x >> 4);
	x |= (x >> 8);
	return (bit_count(x >> 1));
      }
    };
    
    template <>
    struct struct_floor_log2<4>
    {
      uint32_t operator()(uint32_t x) const
      {
	x |= (x >> 1);
        x |= (x >> 2);
        x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
	return (bit_count(x >> 1));
      }
    };
    
    template <>
    struct struct_floor_log2<8>
    {
      uint64_t operator()(uint64_t x) const
      {
	x |= (x >> 1);
        x |= (x >> 2);
        x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);
	x |= (x >> 32);
	return (bit_count(x >> 1));
      }
    };
    
    template <uint64_t X>
    struct static_floor_log2
    {
    private:
      static const uint64_t X1 = uint64_t(X)  | (uint64_t(X)  >>  1);
      static const uint64_t X2 = uint64_t(X1) | (uint64_t(X1) >>  2);
      static const uint64_t X3 = uint64_t(X2) | (uint64_t(X2) >>  4);
      static const uint64_t X4 = uint64_t(X3) | (uint64_t(X3) >>  8);
      static const uint64_t X5 = uint64_t(X4) | (uint64_t(X4) >> 16);
      static const uint64_t X6 = uint64_t(X5) | (uint64_t(X5) >> 32);
      static const uint64_t X7 = X6 >> 1;
    public:
      static const uint64_t result = static_bit_count<X7>::result;
      static const uint64_t value = static_bit_count<X7>::result;
    };
    
    template <typename Tp>
    inline
    uint64_t floor_log2(Tp x)
    {
      static struct_floor_log2<sizeof(Tp)> __func;
      return __func(x);
    }
    
  };
};

#endif
