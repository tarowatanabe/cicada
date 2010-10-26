// -*- mode: c++ -*-

#ifndef __UTILF__BIT_VECTOR__HPP__
#define __UTILF__BIT_VECTOR__HPP__ 1

#include <stdint.h>

#include <iostream>
#include <stdexcept>

#include <utils/bithack.hpp>
#include <utils/hashmurmur.hpp>

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
	0x00000001, 0x00000002, 0x00000004, 0x00000008,
	0x00000010, 0x00000020, 0x00000040, 0x00000080,
	0x00000100, 0x00000200, 0x00000400, 0x00000800,
	0x00001000, 0x00002000, 0x00004000, 0x00008000,
	0x00010000, 0x00020000, 0x00040000, 0x00080000,
	0x00100000, 0x00200000, 0x00400000, 0x00800000,
	0x01000000, 0x02000000, 0x04000000, 0x08000000,
	0x10000000, 0x20000000, 0x40000000, 0x80000000,
      };
      return __mask_blocks;
    }
    
    const block_type* __masks_reverse() const 
    {
      static block_type __mask_blocks[__bitblock_bit_size] = {
	~0x00000001, ~0x00000002, ~0x00000004, ~0x00000008,
	~0x00000010, ~0x00000020, ~0x00000040, ~0x00000080,
	~0x00000100, ~0x00000200, ~0x00000400, ~0x00000800,
	~0x00001000, ~0x00002000, ~0x00004000, ~0x00008000,
	~0x00010000, ~0x00020000, ~0x00040000, ~0x00080000,
	~0x00100000, ~0x00200000, ~0x00400000, ~0x00800000,
	~0x01000000, ~0x02000000, ~0x04000000, ~0x08000000,
	~0x10000000, ~0x20000000, ~0x40000000, ~0x80000000,
      };
      return __mask_blocks;
    }
    
    const block_type* __masks_rank() const 
    {
      static block_type __mask_blocks[__bitblock_bit_size] = {
	0x00000001, 0x00000003, 0x00000007, 0x0000000f,
	0x0000001f, 0x0000003f, 0x0000007f, 0x000000ff,
	0x000001ff, 0x000003ff, 0x000007ff, 0x00000fff,
	0x00001fff, 0x00003fff, 0x00007fff, 0x0000ffff,
	0x0001ffff, 0x0003ffff, 0x0007ffff, 0x000fffff,
	0x001fffff, 0x003fffff, 0x007fffff, 0x00ffffff,
	0x01ffffff, 0x03ffffff, 0x07ffffff, 0x0fffffff,
	0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff,
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
      std::copy((const size_type*) x.begin(), (const size_type*) x.end(), (size_type*) begin());
    }
    
    void clear() 
    {
      std::fill((size_type*) begin(), (size_type*) end(), size_type(0));
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
	return x | y;
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
      __assign_operator_aux((const size_type*) x.begin(), (const size_type*) x.end(), (size_type*) begin(), __or_op<size_type>());
      return *this;
    }
    
    bit_vector& operator&=(const bit_vector& x)
    {
      __assign_operator_aux((const size_type*) x.begin(), (const size_type*) x.end(), (size_type*) begin(), __and_op<size_type>());
      return *this;
    }
    
    bit_vector& operator^=(const bit_vector& x)
    {
      __assign_operator_aux((const size_type*) x.begin(), (const size_type*) x.end(), (size_type*) begin(), __xor_op<size_type>());
      
      return *this;
    }
    
  public:
    size_type count() const
    {
      size_type sum = 0;
      const size_type* biter_end = (const size_type*) end();
      for (const size_type* biter = (const size_type*) begin(); biter != biter_end; ++ biter)
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
    
  private:
    block_type __bitblock[__bitblock_size];
  };
  
  template <size_t _N>
  inline
  std::ostream& operator<<(std::ostream& os, const bit_vector<_N>& x)
  {
    typedef typename bit_vector<_N>::size_type size_type;
    
    for (size_type i = 0; i < x.size(); ++ i) {
      const char mask = (x.test(i) - 1);
      os << char((~mask & '1') | (mask & '0'));
    }

    return os;
  }
  
  template <size_t _N>
  inline
  size_t hash_value(bit_vector<_N> const& x)
  {
    return utils::hashmurmur<size_t>()(x.begin(), x.end(), 0);
  }

  
  template <size_t _N>
  inline
  bool operator==(const bit_vector<_N>& x, const bit_vector<_N>& y)
  {
    return std::equal(x.begin(), x.end(), y.begin());
  }
  
  template <size_t _N>
  inline
  bool operator!=(const bit_vector<_N>& x, const bit_vector<_N>& y)
  {
    return ! (x == y);
  }
  
  template <size_t _N>
  inline
  bool operator<(const bit_vector<_N>& x, const bit_vector<_N>& y)
  {
    return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
  }
  
  template <size_t _N>
  inline
  bool operator>(const bit_vector<_N>& x, const bit_vector<_N>& y)
  {
    return y < x;
  }
  
  template <size_t _N>
  inline
  bool operator<=(const bit_vector<_N>& x, const bit_vector<_N>& y)
  {
    return ! (y < x);
  }
  
  template <size_t _N>
  inline
  bool operator>=(const bit_vector<_N>& x, const bit_vector<_N>& y)
  {
    return ! (x < y);
  }
  

};

#endif
