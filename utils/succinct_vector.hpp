// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SUCCINCT_VECTOR__HPP__
#define __UTILS__SUCCINCT_VECTOR__HPP__ 1

#include <stdint.h>

#include <sstream>
#include <iostream>
#include <vector>
#include <stdexcept>

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/thread.hpp>

#include <utils/atomicop.hpp>
#include <utils/bithack.hpp>
#include <utils/repository.hpp>
#include <utils/map_file.hpp>
#include <utils/array_power2.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/filesystem.hpp>

namespace utils
{
  struct __succinct_vector_base
  {
    typedef boost::filesystem::path path_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef uint32_t block_type;
    typedef uint32_t rank_high_type;
    typedef uint8_t  rank_low_type;
    typedef uint8_t  byte_type;
    
    static const size_type num_block_rank_high = 8;
    static const size_type num_block_rank_low = 1;

    static const size_type block_size     = sizeof(block_type) * 8;           // 32 bit
    static const size_type rank_high_size = block_size * num_block_rank_high; // 256 bit
    static const size_type rank_low_size  = block_size * num_block_rank_low;  // 32 bit
    
    static const size_type shift_block = 5;
    static const size_type shift_block_rank_high = 8;
    static const size_type shift_block_rank_low = 5;

    static const block_type mask_block = 0x1f;
    
    const byte_type* masks_select1() const
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

    const byte_type* masks_select0() const
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
    
    const block_type* masks_rank() const 
    {
      static block_type __mask_blocks[block_size] = {
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

    const block_type* masks() const 
    {
      static block_type __mask_blocks[block_size] = {
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
    
    const block_type* masks_reverse() const 
    {
      static block_type __mask_blocks[block_size] = {
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
    
  private:
    struct rank_block_type
    {
      rank_block_type(const block_type& __block_value) : block_value(__block_value) {}
      
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

    template <typename Data, size_t ReverseScale>
    struct __lower_bound_linear
    {
      static inline
      size_type result(size_type first, size_type last, const Data& data, const size_type& value)
      {
	size_type offset = first;
	for (/**/; first != last && (first - offset + 1) * ReverseScale - data[first] < value; ++ first);
	return first;
      }
    };
    
    template <typename Data>
    struct __lower_bound_linear<Data, 0>
    {
      static inline
      size_type result(size_type first, size_type last, const Data& data, const size_type& value)
      {
	for (/**/; first != last && data[first] < value; ++ first);
	return first;
      }
    };

    template <typename Data, size_t ReverseScale>
    struct __lower_bound
    {
      static inline
      size_type result(size_type first, size_type last, const Data& data, const size_type& value)
      {
	size_type offset = first;
	size_type length = last - first;
	if (length <= 64) {
	  for (/**/; first != last && (first - offset + 1) * ReverseScale - data[first] < value; ++ first);
	  return first;
	} else {
	  while (length > 0) {
	    const size_type half = length >> 1;
	    const size_type middle = first + half;
	    if ((middle - offset + 1) * ReverseScale - data[middle] < value) {
	      first = middle + 1;
	      length = length - half - 1;
	    } else
	      length = half;
	  }
	  return first;
	}
      }
    };
    
    template <typename Data>
    struct __lower_bound<Data, 0>
    {
      static inline
      size_type result(size_type first, size_type last, const Data& data, const size_type& value)
      {
	size_type length = last - first;
	if (length <= 64) {
	  for (/**/; first != last && data[first] < value; ++ first);
	  return first;
	} else {
	  while (length > 0) {
	    const size_type half = length >> 1;
	    const size_type middle = first + half;
	    if (data[middle] < value) {
	      first = middle + 1;
	      length = length - half - 1;
	    } else
	      length = half;
	  }
	  return first;
	}
      }
    };
    
  protected:
    
    template <typename Block, typename RankHigh, typename RankLow>
    size_type select1(const Block& block, const RankHigh& rank_high, const RankLow& rank_low, size_type x) const
    {
      const size_type pos_high = __lower_bound<RankHigh, 0>::result(0, rank_high.size(), rank_high, x);
      if (pos_high == rank_high.size()) return size_type(-1);
      size_type num_bits_remain = x - (pos_high == 0 ? size_type(0) : rank_high[pos_high - 1]);
      
      const size_type pos_low_first = pos_high * 7;
      const size_type pos_low_last = std::min(pos_low_first + 7, rank_low.size());
      const size_type pos_low = __lower_bound_linear<RankLow, 0>::result(pos_low_first, pos_low_last, rank_low, num_bits_remain);
      const size_type pos_low_diff = pos_low - pos_low_first;
      num_bits_remain -= (pos_low_diff == 0 ? size_type(0) : rank_low[pos_low - 1]);
      
      // now, pos_low is the position of block...
      const size_type pos_block = pos_high * num_block_rank_high + pos_low_diff;
      const rank_block_type rank_block(block[pos_block]);
      const size_type pos_byte = __lower_bound_linear<rank_block_type, 0>::result(0, rank_block.size() - 1, rank_block, num_bits_remain);
      num_bits_remain -= (pos_byte == 0 ? size_type(0) : rank_block[pos_byte - 1]);
      
      size_type pos_bit = (pos_block << shift_block) + (pos_byte << 3);
      const byte_type byte_value = rank_block(pos_byte);
      const size_type byte_value_count_lower = bithack::bit_count(byte_value & 0x0f);
      
      pos_bit += (byte_value_count_lower < num_bits_remain
		  ? 4 + masks_select1()[4 * ((byte_value >> 4) & 0x0f) + (num_bits_remain - byte_value_count_lower - 1)]
		  : masks_select1()[4 * (byte_value & 0x0f) + (num_bits_remain - 1)]);
      
      return pos_bit;
    }
    
    template <typename Block, typename RankHigh, typename RankLow>
    size_type select0(const Block& block, const RankHigh& rank_high, const RankLow& rank_low, size_type x) const
    {
      const size_type pos_high = __lower_bound<RankHigh, rank_high_size>::result(0, rank_high.size(), rank_high, x);
      if (pos_high == rank_high.size()) return size_type(-1);
      size_type num_bits_remain = x - (pos_high == 0 ? size_type(0) : (pos_high << shift_block_rank_high) - rank_high[pos_high - 1]);
      
      const size_type pos_low_first = pos_high * 7;
      const size_type pos_low_last = std::min(pos_low_first + 7, rank_low.size());
      const size_type pos_low = __lower_bound_linear<RankLow, 32>::result(pos_low_first, pos_low_last, rank_low, num_bits_remain);
      const size_type pos_low_diff = pos_low - pos_low_first;
      num_bits_remain -= (pos_low_diff == 0 ? size_type(0) : (pos_low_diff << shift_block_rank_low) - rank_low[pos_low - 1]);
      
      // now, pos_low is the position of block...
      const size_type pos_block = pos_high * num_block_rank_high + pos_low_diff;
      const rank_block_type rank_block(block[pos_block]);
      const size_type pos_byte = __lower_bound_linear<rank_block_type, 8>::result(0, rank_block.size() - 1, rank_block, num_bits_remain);
      num_bits_remain -= (pos_byte == 0 ? size_type(0) : (pos_byte << 3) - rank_block[pos_byte - 1]);
      
      size_type pos_bit = (pos_block << shift_block) + (pos_byte << 3);
      const byte_type byte_value = rank_block(pos_byte);
      const size_type byte_value_count_lower = 4 - bithack::bit_count(byte_value & 0x0f);
      
      pos_bit += (byte_value_count_lower < num_bits_remain
		  ? 4 + masks_select0()[4 * ((byte_value >> 4) & 0x0f) + (num_bits_remain - byte_value_count_lower - 1)]
		  : masks_select0()[4 * (byte_value & 0x0f) + (num_bits_remain - 1)]);
      
      return pos_bit;
    }
    
    template <typename Block, typename RankHigh, typename RankLow>
    size_type rank1(const Block& block, const RankHigh& rank_high, const RankLow& rank_low, size_type pos) const
    {
      const size_type pos_rank_high = pos >> shift_block_rank_high;
      const size_type pos_rank_low  = pos >> shift_block_rank_low;
      const size_type pos_rank_low_offset = pos_rank_low & (num_block_rank_high - 1);
      
      // we need to tweak this!
      return ((pos_rank_high == 0 ? size_type(0) : rank_high[pos_rank_high - 1])
	      + (pos_rank_low_offset == 0 ? size_type(0) : rank_low[pos_rank_high * 7 + pos_rank_low_offset - 1])
	      + bithack::bit_count(block[pos_rank_low] & masks_rank()[pos & mask_block]));
    }    
  };

  template <typename _Alloc=std::allocator<uint32_t> >
  class succinct_vector_mapped : protected __succinct_vector_base
  {
  private:
    typedef typename _Alloc::template rebind<block_type>::other      block_allocator_type;
    typedef typename _Alloc::template rebind<rank_high_type>::other  rank_high_allocator_type;
    typedef typename _Alloc::template rebind<rank_low_type>::other   rank_low_allocator_type;
    
    typedef utils::map_file<block_type, block_allocator_type>  bit_block_type;
    typedef utils::map_file<rank_high_type, rank_high_allocator_type> bit_rank_high_type;
    typedef utils::map_file<rank_low_type, rank_low_allocator_type>   bit_rank_low_type;
    
    typedef __succinct_vector_base base_type;
    typedef succinct_vector_mapped<_Alloc> self_type;
    
  public:
    typedef base_type::size_type       size_type;
    typedef base_type::difference_type difference_type;

  private:
    struct Cache
    {
      typedef int64_t value_type;
      
      Cache() : value(value_type(-1)) {}
      
      volatile value_type value;
    };
    typedef Cache cache_type;
    
    typedef typename _Alloc::template rebind<cache_type>::other cache_allocator_type;
    typedef std::vector<cache_type, cache_allocator_type > cache_set_type;

  public:
    succinct_vector_mapped()
      : __size(0), __blocks(), __rank_high(), __rank_low() {}
    succinct_vector_mapped(const path_type& path)
      : __size(0), __blocks(), __rank_high(), __rank_low() { open(path); }
    
    bool operator[](size_type pos) const { return test(pos); }
    bool test(size_type pos) const
    {
      return __blocks[pos >> shift_block] & masks()[pos & mask_block];
    }
    size_type size() const { return __size; }
    bool empty() const { return __size == 0; }
    bool is_open() const { return __blocks.is_open(); }
    
    uint64_t size_bytes() const { return __blocks.size_bytes() + __rank_high.size_bytes() + __rank_low.size_bytes(); }
    uint64_t size_compressed() const { return __blocks.size_compressed() + __rank_high.size_compressed() + __rank_low.size_compressed(); }
    uint64_t size_cache() const { return __cache_select0.size() * sizeof(cache_type) + __cache_select1.size() * sizeof(cache_type); }
    
    void clear()
    {
      __size = 0;
      __blocks.clear();
      __rank_high.clear();
      __rank_low.clear();
      
      __cache_select0.clear();
      __cache_select1.clear();
    }
    void close() { clear(); }
    
    size_type select(size_type pos, bool bit) const
    {
      if (__rank_high.empty() || __rank_low.empty()) 
	throw std::runtime_error("no ranks...");

      cache_set_type& caches = const_cast<cache_set_type&>(bit ? __cache_select1 : __cache_select0);
      const uint64_t mask_cache = caches.size() - 1;
      const uint64_t mask_pos    = (~uint64_t(bit - 1) & __select1_mask_pos)    | (uint64_t(bit - 1) & __select0_mask_pos);
      const uint64_t mask_select = (~uint64_t(bit - 1) & __select1_mask_select) | (uint64_t(bit - 1) & __select0_mask_select);
      
      cache_type cache;
      cache_type cache_new;
      
      cache.value = utils::atomicop::fetch_and_add(caches[pos & mask_cache].value, int64_t(0));
      
      uint64_t __pos    = (cache.value & mask_pos);
      uint64_t __select = (cache.value & mask_select);
      
      if (__pos == ((uint64_t(pos) << 32) & mask_pos))
	return __select;
      
      __select = (bit 
		  ? base_type::select1(__blocks, __rank_high, __rank_low, pos)
		  : base_type::select0(__blocks, __rank_high, __rank_low, pos));
      
      cache_new.value = ((uint64_t(pos) << 32) & mask_pos) | (uint64_t(__select) & mask_select);
      
      utils::atomicop::compare_and_swap(caches[pos & mask_cache].value, cache.value, cache_new.value);
      
      return __select;
    }
    
    size_type rank(size_type pos, bool bit) const
    {
      if (__rank_high.empty() || __rank_low.empty()) 
	throw std::runtime_error("no ranks...");

      const size_type rank1_value = base_type::rank1(__blocks, __rank_high, __rank_low, pos);
      return (~size_type(bit - 1) & rank1_value) | (size_type(bit - 1) & (pos + 1 - rank1_value));
    }
    
    path_type path() const { return __blocks.path().parent_path(); }

    static bool exists(const path_type& path) 
    {
      if (! utils::repository::exists(path)) return false;
      if (! bit_block_type::exists(path / "bits")) return false;
      if (! bit_rank_high_type::exists(path / "rank-high")) return false;
      if (! bit_rank_low_type::exists(path / "rank-low")) return false;
      return true;
    }

    void read(const path_type& path) { open(path); }
    void open(const path_type& path)
    {
      typedef utils::repository repository_type;
      
      close();
      
      repository_type repository(path, repository_type::read);
      __blocks.open(repository.path("bits"));
      __rank_high.open(repository.path("rank-high"));
      __rank_low.open(repository.path("rank-low"));
      
      repository_type::const_iterator iter = repository.find("size");
      if (iter == repository.end())
	throw std::runtime_error("no size...");
      __size = atoll(iter->second.c_str());

      repository_type::const_iterator titer = repository.find("type");
      if (titer == repository.end())
	throw std::runtime_error("no type...");
      if (titer->second != "succinct-vector")
	throw std::runtime_error("not a succinct vector...");
      
      // setup cache...
      const size_type select0_size = std::max(size_type(utils::bithack::next_largest_power2((size() - __rank_high.back()) >> 10)),
					      size_type(1024 * 32));
      const size_type select1_size = std::max(size_type(utils::bithack::next_largest_power2(__rank_high.back() >> 10)),
					      size_type(1024 * 32));
      
      __cache_select0.reserve(select0_size);
      __cache_select0.resize(select0_size, cache_type());
      
      __cache_select1.reserve(select1_size);
      __cache_select1.resize(select1_size, cache_type());
      
      __select0_mask_pos = (~uint64_t(__cache_select0.size() - 1)) << 32;
      __select1_mask_pos = (~uint64_t(__cache_select1.size() - 1)) << 32;
      __select0_mask_select = ~__select0_mask_pos;
      __select1_mask_select = ~__select1_mask_pos;
    }
    
    void write(const path_type& file) const
    {
      if (path() == file) return;
      
      // remove first...
      if (boost::filesystem::exists(file) && ! boost::filesystem::is_directory(file))
	boost::filesystem::remove_all(file);
      
      // create directory
      if (! boost::filesystem::exists(file))
	boost::filesystem::create_directories(file);
      
      // remove all the files...
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(file); iter != iter_end; ++ iter)
	boost::filesystem::remove_all(*iter);
      
      // copy all...
      for (boost::filesystem::directory_iterator iter(path()); iter != iter_end; ++ iter)
	utils::filesystem::copy_files(*iter, file);
    }
    
  public:
    size_type          __size;
    bit_block_type     __blocks;
    bit_rank_high_type __rank_high;
    bit_rank_low_type  __rank_low;
    
  private:
    cache_set_type __cache_select0;
    cache_set_type __cache_select1;

    uint64_t __select0_mask_pos;
    uint64_t __select0_mask_select;
    uint64_t __select1_mask_pos;
    uint64_t __select1_mask_select;
  };
  
  template <typename _Alloc=std::allocator<uint32_t> >
  class succinct_vector : protected __succinct_vector_base
  {
  private:
    typedef typename _Alloc::template rebind<block_type>::other      block_allocator_type;
    typedef typename _Alloc::template rebind<rank_high_type>::other  rank_high_allocator_type;
    typedef typename _Alloc::template rebind<rank_low_type>::other   rank_low_allocator_type;
    
    typedef std::vector<block_type, block_allocator_type>         bit_block_type;
    typedef std::vector<rank_high_type, rank_high_allocator_type> bit_rank_high_type;
    typedef std::vector<rank_low_type, rank_low_allocator_type>   bit_rank_low_type;

    typedef __succinct_vector_base base_type;

  public:
    typedef base_type::size_type       size_type;
    typedef base_type::difference_type difference_type;
    
  public:
    succinct_vector()
      : __size(0),
	__blocks(),
	__rank_high(),
	__rank_low() {}
    
    template <typename __Alloc>
    succinct_vector(const succinct_vector_mapped<__Alloc>& x)
      : __size(x.__size),
	__blocks(x.__blocks.begin(), x.__blocks.end()),
	__rank_high(x.__rank_high.begin(), x.__rank_high.end()),
	__rank_low(x.__rank_low.begin(), x.__rank_low.end()) {}
    
    template <typename __Alloc>
    succinct_vector& operator=(const succinct_vector_mapped<__Alloc>& x)
    {
      __size = x.__size;
      __blocks.assign(x.__blocks.begin(), x.__blocks.end());
      __rank_high.assign(x.__rank_high.begin(), x.__rank_high.end());
      __rank_low.assign(x.__rank_low.begin(), x.__rank_low.end());
      return *this;
    }
    
    void set(size_type pos, bool bit=true)
    {
      __rank_high.clear();
      __rank_low.clear();

      const size_type block_pos = pos >> shift_block;
      const size_type mask_pos  = pos & mask_block;
      
      __size = bithack::max(pos + 1, __size);
      if (block_pos >= __blocks.size())
	__blocks.resize(block_pos + 1, block_type());
      
      __blocks[block_pos] = (__blocks[block_pos] & masks_reverse()[mask_pos]) | (-bit & masks()[mask_pos]);
    }

    bool operator[](size_type pos) const { return test(pos); }
    bool test(size_type pos) const
    {
      return __blocks[pos >> shift_block] & masks()[pos & mask_block];
    }
    size_type size() const { return __size; }
    bool empty() const { return __size == 0; }

    uint64_t size_bytes() const
    { 
      return (__blocks.size() * sizeof(block_type)
	      + __rank_high.size() * sizeof(rank_high_type)
	      + __rank_low.size() * sizeof(rank_low_type));
    }
    uint64_t size_compressed() const { return size_bytes(); }
    uint64_t size_cache() const { return 0; }
    
    void clear()
    {
      __size = 0;
      __blocks.clear();
      __rank_high.clear();
      __rank_low.clear();
    }
    
    size_type select(size_type pos, bool bit) const
    {
      if (__rank_high.empty() || __rank_low.empty()) 
	throw std::runtime_error("no ranks...");
      
      return (bit 
	      ? base_type::select1(__blocks, __rank_high, __rank_low, pos)
	      : base_type::select0(__blocks, __rank_high, __rank_low, pos));
    }
    
    size_type rank(size_type pos, bool bit) const
    {
      if (__rank_high.empty() || __rank_low.empty()) 
	throw std::runtime_error("no ranks...");
      
      const size_type rank1_value = base_type::rank1(__blocks, __rank_high, __rank_low, pos);
      return (~size_type(bit - 1) & rank1_value) | (size_type(bit - 1) & (pos + 1 - rank1_value));
    }

    
    void build()
    {
      __rank_high.clear();
      __rank_low.clear();
      
      rank_high_type sum = 0;
      rank_high_type sum_low = 0;
      
      const size_type mask_dump = num_block_rank_high - 1;
      
      // we emit hight block every num_block_rank_high
      int i = 0;
      typename bit_block_type::const_iterator biter_end = __blocks.end();
      for (typename bit_block_type::const_iterator biter = __blocks.begin(); biter != biter_end; ++ biter, ++ i) {
	const rank_high_type bitcount = bithack::bit_count(*biter);
	
	sum_low += bitcount;
	sum += bitcount;
	
	if ((i & mask_dump) == mask_dump) {
	  __rank_high.push_back(sum);
	  sum_low = 0;
	} else
	  __rank_low.push_back(sum_low);
      }
      
      if ((i & mask_dump) != 0 || sum == 0)
	__rank_high.push_back(sum);
    }

    void write(const path_type& path) const
    {
      typedef utils::repository repository_type;
      
      if (__rank_high.empty())
	const_cast<succinct_vector&>(*this).build();

      repository_type repository(path, repository_type::write);
      dump_file(repository.path("bits"), __blocks);
      dump_file(repository.path("rank-high"), __rank_high);
      dump_file(repository.path("rank-low"), __rank_low);
      
      std::ostringstream stream;
      stream << __size;
      repository["size"] = stream.str();
      repository["type"] = "succinct-vector";
    }

  private:
    template <typename _Path, typename _Data>
    inline
    void dump_file(const _Path& file, const _Data& data) const
    {
      boost::iostreams::filtering_ostream os;
      os.push(boost::iostreams::file_sink(file.native_file_string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
      
      const int64_t file_size = sizeof(typename _Data::value_type) * data.size();
      for (int64_t offset = 0; offset < file_size; offset += 1024 * 1024)
	if (! os.write(((char*) &(*data.begin())) + offset, std::min(int64_t(1024 * 1024), file_size - offset)))
	  throw std::runtime_error("succinct_vector write():");
    }
    
    
  private:
    size_type          __size;
    bit_block_type     __blocks;
    bit_rank_high_type __rank_high;
    bit_rank_low_type  __rank_low;
  };
  
};


#endif
