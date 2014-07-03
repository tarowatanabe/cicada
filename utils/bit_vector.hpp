// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILF__BIT_VECTOR__HPP__
#define __UTILF__BIT_VECTOR__HPP__ 1

#include <stdint.h>

#include <iostream>
#include <stdexcept>

#include <utils/bithack.hpp>
#include <utils/hashmurmur3.hpp>

namespace utils
{

  struct __bit_vector_base
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef uint32_t  block_type;
    typedef uint8_t   byte_type;
    
  public:
    static const size_type __bitblock_byte_size = sizeof(block_type);
    static const size_type __bitblock_bit_size  = __bitblock_byte_size * 8;
    static const size_type __bitblock_mask      = __bitblock_bit_size - 1;
    static const size_type __bitblock_shift     = utils::bithack::static_bit_count<__bitblock_mask>::result;

    const block_type* __masks() const 
    {
      static block_type __mask_blocks[__bitblock_bit_size] = {
	0x00000001u, 0x00000002u, 0x00000004u, 0x00000008u,
	0x00000010u, 0x00000020u, 0x00000040u, 0x00000080u,
	0x00000100u, 0x00000200u, 0x00000400u, 0x00000800u,
	0x00001000u, 0x00002000u, 0x00004000u, 0x00008000u,
	0x00010000u, 0x00020000u, 0x00040000u, 0x00080000u,
	0x00100000u, 0x00200000u, 0x00400000u, 0x00800000u,
	0x01000000u, 0x02000000u, 0x04000000u, 0x08000000u,
	0x10000000u, 0x20000000u, 0x40000000u, 0x80000000u,
      };
      return __mask_blocks;
    }
    
    const block_type* __masks_reverse() const 
    {
      static block_type __mask_blocks[__bitblock_bit_size] = {
	~0x00000001u, ~0x00000002u, ~0x00000004u, ~0x00000008u,
	~0x00000010u, ~0x00000020u, ~0x00000040u, ~0x00000080u,
	~0x00000100u, ~0x00000200u, ~0x00000400u, ~0x00000800u,
	~0x00001000u, ~0x00002000u, ~0x00004000u, ~0x00008000u,
	~0x00010000u, ~0x00020000u, ~0x00040000u, ~0x00080000u,
	~0x00100000u, ~0x00200000u, ~0x00400000u, ~0x00800000u,
	~0x01000000u, ~0x02000000u, ~0x04000000u, ~0x08000000u,
	~0x10000000u, ~0x20000000u, ~0x40000000u, ~0x80000000u,
      };
      return __mask_blocks;
    }
    
    const block_type* __masks_rank() const 
    {
      static block_type __mask_blocks[__bitblock_bit_size] = {
	0x00000001u, 0x00000003u, 0x00000007u, 0x0000000fu,
	0x0000001fu, 0x0000003fu, 0x0000007fu, 0x000000ffu,
	0x000001ffu, 0x000003ffu, 0x000007ffu, 0x00000fffu,
	0x00001fffu, 0x00003fffu, 0x00007fffu, 0x0000ffffu,
	0x0001ffffu, 0x0003ffffu, 0x0007ffffu, 0x000fffffu,
	0x001fffffu, 0x003fffffu, 0x007fffffu, 0x00ffffffu,
	0x01ffffffu, 0x03ffffffu, 0x07ffffffu, 0x0fffffffu,
	0x1fffffffu, 0x3fffffffu, 0x7fffffffu, 0xffffffffu,
      };
      return __mask_blocks;
    }

    const byte_type* __masks_select1() const
    {
      static byte_type __mask_blocks[16 * 4] = {
	/* 3210 */
	/* 0000 */ 0, 0, 0, 0,
	/* 0001 */ 0, 0, 0, 0,
	/* 0010 */ 1, 0, 0, 0,
	/* 0011 */ 0, 1, 0, 0,
	/* 0100 */ 2, 0, 0, 0,
	/* 0101 */ 0, 2, 0, 0,
	/* 0110 */ 1, 2, 0, 0,
	/* 0111 */ 0, 1, 2, 0,
	/* 1000 */ 3, 0, 0, 0,
	/* 1001 */ 0, 3, 0, 0,
	/* 1010 */ 1, 3, 0, 0,
	/* 1011 */ 0, 1, 3, 0,
	/* 1100 */ 2, 3, 0, 0,
	/* 1101 */ 0, 2, 3, 0,
	/* 1110 */ 1, 2, 3, 0,
	/* 1111 */ 0, 1, 2, 3,
      };
      return __mask_blocks;
    }

    const byte_type* __masks_select0() const
    {
      static byte_type __mask_blocks[16 * 4] = {
	/* 3210 */
	/* 0000 */ 0, 1, 2, 3,
	/* 0001 */ 1, 2, 3, 0,
	/* 0010 */ 0, 2, 3, 0,
	/* 0011 */ 2, 3, 0, 0,
	/* 0100 */ 0, 1, 3, 0,
	/* 0101 */ 1, 3, 0, 0,
	/* 0110 */ 0, 3, 0, 0,
	/* 0111 */ 3, 0, 0, 0,
	/* 1000 */ 0, 1, 2, 0,
	/* 1001 */ 1, 2, 0, 0,
	/* 1010 */ 0, 2, 0, 0,
	/* 1011 */ 2, 0, 0, 0,
	/* 1100 */ 0, 1, 0, 0,
	/* 1101 */ 1, 0, 0, 0,
	/* 1110 */ 0, 0, 0, 0,
	/* 1111 */ 0, 0, 0, 0,
      };
      return __mask_blocks;
    }
  };

  template <typename BitVector>
  struct bit_vector_expr_or : public __bit_vector_base
  {
    bit_vector_expr_or(const BitVector& x, const BitVector& y)
      : bv1(x), bv2(y) {}

    size_type count() const
    {
      size_type sum = 0;

      const block_type* iter1 = (const block_type*) bv1.begin();
      const block_type* iter1_end = (const block_type*) bv1.end();

      const block_type* iter2 = (const block_type*) bv2.begin();
      for (/**/; iter1 != iter1_end; ++ iter1, ++ iter2)
	sum += utils::bithack::bit_count(*iter1 | *iter2);
      
      return sum;
    }
    
    const BitVector& bv1;
    const BitVector& bv2;
  };

  
  template <size_t NumBits>
  class bit_vector : public __bit_vector_base
  {
  public:
    static const size_type __bit_multiple_size = 128;
    static const size_type __bit_size          = NumBits;
    static const size_type __bit_capacity      = (__bit_size + (__bit_multiple_size - 1)) & size_type(- __bit_multiple_size);    
    
    static const size_type __bitblock_size      = __bit_capacity >> __bitblock_shift;
    
  public:
    bit_vector() { clear(); }
    bit_vector(const bit_vector& x) { assign(x); }
    bit_vector& operator=(const bit_vector& x)
    {
      assign(x);
      return *this;
    }

  public:
    size_type size() const { return __bit_size; }
    size_type capacity() const { return __bit_capacity; }
    size_type block_size() const { return __bitblock_size; }
    
  public:
    void assign(const bit_vector& x)
    {
      std::copy((const block_type*) x.begin(), (const block_type*) x.end(), (block_type*) begin());
    }
    
    void clear() 
    {
      std::fill((block_type*) begin(), (block_type*) end(), block_type(0));
    }

  public:
    bool operator[](size_type pos) const { return test(pos); }
    
    bool test(size_type pos) const
    {
      return __bitblock[pos >> __bitblock_shift] & __masks()[pos & __bitblock_mask];
    }
    
    void clear(size_type pos) { set(pos, false); }
    
    void set(size_type pos, bool bit=true)
    {
      const size_type pos_block = pos >> __bitblock_shift;
      const size_type pos_mask  = pos & __bitblock_mask;
      
      __bitblock[pos_block] = (__bitblock[pos_block] & __masks_reverse()[pos_mask]) | (-bit & __masks()[pos_mask]);
    }
    
  private:    
    template <typename Tp>
    struct __or_op
    {
      Tp operator()(const Tp& x, const Tp& y) const
      {
	return x | y;
      }
    };
    
    template <typename Tp>
    struct __and_op
    {
      Tp operator()(const Tp& x, const Tp& y) const
      {
	return x & y;
      }
    };
    
    template <typename Tp>
    struct __xor_op
    {
      Tp operator()(const Tp& x, const Tp& y) const
      {
	return x ^ y;
      }
    };
    
    template <typename Iterator1, typename Iterator2, typename Operator>
    void __assign_operator_aux(Iterator1 first1, Iterator1 last1, Iterator2 first2, Operator op)
    {
      for (/**/; first1 != last1; ++ first1, ++ first2)
	*first2 = op(*first2, *first1);
    }
    
  public:
    bit_vector& operator|=(const bit_vector& x)
    {
      __assign_operator_aux((const block_type*) x.begin(), (const block_type*) x.end(), (block_type*) begin(), __or_op<block_type>());
      return *this;
    }
    
    bit_vector& operator&=(const bit_vector& x)
    {
      __assign_operator_aux((const block_type*) x.begin(), (const block_type*) x.end(), (block_type*) begin(), __and_op<block_type>());
      return *this;
    }
    
    bit_vector& operator^=(const bit_vector& x)
    {
      __assign_operator_aux((const block_type*) x.begin(), (const block_type*) x.end(), (block_type*) begin(), __xor_op<block_type>());
      
      return *this;
    }
    
  public:
    size_type count() const
    {
      size_type sum = 0;
      const block_type* biter_end = (const block_type*) end();
      for (const block_type* biter = (const block_type*) begin(); biter != biter_end; ++ biter)
	sum += utils::bithack::bit_count(*biter);
      return sum;
    }
    
    size_type count(size_type bits) const
    {
      return (bits == 0 ? size_type(0) : rank(bits - 1, 1));
    }

    size_type rank(size_type pos, bool bit) const
    {
      const size_type pos_block = pos >> __bitblock_shift;
      const size_type pos_mask  = pos & __bitblock_mask;
      
      // first, computa rank1
      size_type sum = 0;
      for (const block_type* biter = begin(); biter != begin() + pos_block; ++ biter)
	sum += utils::bithack::bit_count(*biter);
      sum += utils::bithack::bit_count(__bitblock[pos_block] & __masks_rank()[pos_mask]);
      
      // then, compute rank1 or rank0 according to bit
      const size_type rank_mask = size_type(bit - 1);
      return (~rank_mask & sum) | (rank_mask & (pos + 1 - sum));
    }


    size_type select(size_type x, bool bit) const
    {
      return (bit ? select1(x) : select0(x));
    }

  private:
    struct rank_block_type
    {
      rank_block_type(const block_type& __bitblock_value) : block_value(__bitblock_value) {}
      
      size_type size() const { return sizeof(block_type); }
      
      byte_type operator()(size_type pos) const
      {
	switch (pos) {
	case 0: return (block_value      ) & 0xff;
	case 1: return (block_value >>  8) & 0xff;
	case 2: return (block_value >> 16) & 0xff;
	case 3: return (block_value >> 24) & 0xff;
	default:
	  throw std::runtime_error("index out of range");
	}
      }
      
      // rank access...
      size_type operator[](size_type pos) const
      {
	switch (pos) {
	case 0: return bithack::bit_count(block_value & 0x000000ff);
	case 1: return bithack::bit_count(block_value & 0x0000ffff);
	case 2: return bithack::bit_count(block_value & 0x00ffffff);
	case 3: return bithack::bit_count(block_value & 0xffffffff);
	default:
	  throw std::runtime_error("index out of range");
	}
      }
      
    private:
      block_type block_value;
    };
    
    struct rank_block_reverse_type
    {
      rank_block_reverse_type(const block_type& __bitblock_value) : block_value(__bitblock_value) {}
      
      size_type size() const { return sizeof(block_type); }
      
      byte_type operator()(size_type pos) const
      {
	switch (pos) {
	case 0: return (block_value      ) & 0xff;
	case 1: return (block_value >>  8) & 0xff;
	case 2: return (block_value >> 16) & 0xff;
	case 3: return (block_value >> 24) & 0xff;
	default:
	  throw std::runtime_error("index out of range");
	}
      }
      
      // rank access...
      size_type operator[](size_type pos) const
      {
	switch (pos) {
	case 0: return  8 - bithack::bit_count(block_value & 0x000000ff);
	case 1: return 16 - bithack::bit_count(block_value & 0x0000ffff);
	case 2: return 24 - bithack::bit_count(block_value & 0x00ffffff);
	case 3: return 32 - bithack::bit_count(block_value & 0xffffffff);
	default:
	  throw std::runtime_error("index out of range");
	}
      }
      
    private:
      block_type block_value;
    };

    size_type select1(size_type x) const
    {
      if (x == 0) return size_type(-1);  // undefined... but who cares select1 for zero?
      
      size_type sum_prev = 0;
      size_type sum = 0;
      const block_type* biter = begin();
      for (/**/; biter != end(); ++ biter) {
	sum_prev = sum;
	sum += utils::bithack::bit_count(*biter);
	
	if (! (sum < x)) break;
      }
      if (biter == end()) return size_type(-1);
      size_type bits_remain = x - sum_prev;
      
      const rank_block_type rank_block(*biter);
      size_type pos_byte = 0;
      for (/**/; pos_byte < 4 && rank_block[pos_byte] < bits_remain; ++ pos_byte);
      bits_remain -= (pos_byte == 0 ? size_type(0) : rank_block[pos_byte - 1]);
      
      size_type pos_bit = ((biter - begin()) << __bitblock_shift) + (pos_byte << 3);
      
      const byte_type byte_value = rank_block(pos_byte);
      const size_type byte_value_count_lower = utils::bithack::bit_count(byte_value & 0x0f);
      
      pos_bit += (byte_value_count_lower < bits_remain
		  ? 4 + __masks_select1()[4 * ((byte_value >> 4) & 0x0f) + (bits_remain - byte_value_count_lower - 1)]
		  : __masks_select1()[4 * (byte_value & 0x0f) + (bits_remain - 1)]);
      
      return pos_bit;
    }
    
    size_type select0(size_type x) const
    {
      if (x == 0) return size_type(-1);  // undefined... but who cares select0 for zero?
      
      size_type sum_prev = 0;
      size_type sum = 0;
      const block_type* biter = begin();
      for (/**/; biter != end(); ++ biter) {
	sum_prev = sum;
	sum += __bitblock_bit_size - utils::bithack::bit_count(*biter);
	
	if (! (sum < x)) break;
      }
      
      if (biter == end()) return size_type(-1);
      size_type bits_remain = x - sum_prev;
      
      const rank_block_reverse_type rank_block(*biter);
      size_type pos_byte = 0;
      for (/**/; pos_byte < 4 && rank_block[pos_byte] < bits_remain; ++ pos_byte);
      bits_remain -= (pos_byte == 0 ? size_type(0) : rank_block[pos_byte - 1]);
      
      size_type pos_bit = ((biter - begin()) << __bitblock_shift) + (pos_byte << 3);
      
      const byte_type byte_value = rank_block(pos_byte);
      const size_type byte_value_count_lower = 4 - bithack::bit_count(byte_value & 0x0f);
      
      pos_bit += (byte_value_count_lower < bits_remain
		  ? 4 + __masks_select0()[4 * ((byte_value >> 4) & 0x0f) + (bits_remain - byte_value_count_lower - 1)]
		  : __masks_select0()[4 * (byte_value & 0x0f) + (bits_remain - 1)]);
      
      return pos_bit;
    }

  public:
    inline const block_type* begin() const { return __bitblock; }
    inline       block_type* begin()       { return __bitblock; }
    
    inline const block_type* end() const { return __bitblock + __bitblock_size; }
    inline       block_type* end()       { return __bitblock + __bitblock_size; }
    
    inline const block_type* end(size_type bits) const { return __bitblock + ((bits + __bitblock_mask) >> __bitblock_shift); }
    inline       block_type* end(size_type bits)       { return __bitblock + ((bits + __bitblock_mask) >> __bitblock_shift); }

  public:
    friend
    bit_vector_expr_or<bit_vector<NumBits> > operator|(const bit_vector<NumBits>& x, const bit_vector<NumBits>& y)
    {
      return bit_vector_expr_or<bit_vector<NumBits> >(x, y);
    }
    
  private:
    block_type __bitblock[__bitblock_size];
  };
  
  template <size_t N>
  inline
  std::ostream& operator<<(std::ostream& os, const bit_vector<N>& x)
  {
    typedef typename bit_vector<N>::size_type size_type;
    
    for (size_type i = 0; i < x.size(); ++ i) {
      const char mask = (x.test(i) - 1);
      os << char((~mask & '1') | (mask & '0'));
    }

    return os;
  }
  
  template <size_t N>
  inline
  size_t hash_value(bit_vector<N> const& x)
  {
    return utils::hashmurmur3<size_t>()(x.begin(), x.end(), 0);
  }

  
  template <size_t N>
  inline
  bool operator==(const bit_vector<N>& x, const bit_vector<N>& y)
  {
    return std::equal(x.begin(), x.end(), y.begin());
  }
  
  template <size_t N>
  inline
  bool operator!=(const bit_vector<N>& x, const bit_vector<N>& y)
  {
    return ! (x == y);
  }
  
  template <size_t N>
  inline
  bool operator<(const bit_vector<N>& x, const bit_vector<N>& y)
  {
    return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
  }
  
  template <size_t N>
  inline
  bool operator>(const bit_vector<N>& x, const bit_vector<N>& y)
  {
    return y < x;
  }
  
  template <size_t N>
  inline
  bool operator<=(const bit_vector<N>& x, const bit_vector<N>& y)
  {
    return ! (y < x);
  }
  
  template <size_t N>
  inline
  bool operator>=(const bit_vector<N>& x, const bit_vector<N>& y)
  {
    return ! (x < y);
  }
  

};

#endif
