// -*- mode: c++ -*-

#ifndef __UTILS__BITPACK__HPP__
#define __UTILS__BITPACK__HPP__ 1

//
// bit-packer using template metaprogramming...
// 
// we pack values into bit-vector of multiple of bit-size-width size vector
//

#include <stdint.h>

#include <algorithm>
#include <stdexcept>

#include <utils/bithack.hpp>

namespace utils
{
  namespace bitpack
  {
    void __pack(const uint8_t*  source, uint8_t*  destination, size_t bit_size);
    void __pack(const uint16_t* source, uint16_t* destination, size_t bit_size);
    void __pack(const uint32_t* source, uint32_t* destination, size_t bit_size);
    void __pack(const uint64_t* source, uint64_t* destination, size_t bit_size);
    
    void __unpack(const uint8_t*  source, uint8_t*  destination, size_t bit_size);
    void __unpack(const uint16_t* source, uint16_t* destination, size_t bit_size);
    void __unpack(const uint32_t* source, uint32_t* destination, size_t bit_size);
    void __unpack(const uint64_t* source, uint64_t* destination, size_t bit_size);

    template <size_t ByteSize>
    struct __unpack_value {};
    template <>
    struct __unpack_value<1> { typedef uint8_t value_type; };
    template <>
    struct __unpack_value<2> { typedef uint16_t value_type; };
    template <>
    struct __unpack_value<4> { typedef uint32_t value_type; };
    template <>
    struct __unpack_value<8> { typedef uint64_t value_type; };
    
    template <typename Tp>
    Tp  __unpack(const Tp*  source, size_t pos, size_t bit_size)
    {
      const size_t unit_size = sizeof(Tp) * 8;
      const size_t bits_total = unit_size * bit_size;
      const Tp     mask = (Tp(1) << bit_size) - 1;
      
      const size_t shift_size_total = bits_total - bit_size * (pos + 1);
      const size_t code_pos = bit_size - 1 - (shift_size_total >> utils::bithack::static_floor_log2<sizeof(Tp)*8>::result);
      const size_t shift_amount = shift_size_total & (unit_size - 1);
      
      typedef typename __unpack_value<sizeof(Tp)>::value_type uvalue_type;
      
      if (unit_size - shift_amount < bit_size)
	return (((uvalue_type(source[code_pos]) >> shift_amount) & mask) 
		| ((source[code_pos - 1] << (unit_size - shift_amount)) & mask));
      else
	return (uvalue_type(source[code_pos]) >> shift_amount) & mask;
    }

    template <typename Tp>
    inline
    void pack(const Tp* source, Tp* destination, size_t size, size_t bit_size)
    {
      typedef typename __unpack_value<sizeof(Tp)>::value_type uvalue_type;

      if (bit_size == sizeof(Tp) * 8)
	std::copy(source, source + size, destination);
      else
	for (size_t i = 0; i < size; i += sizeof(Tp) * 8, source += sizeof(Tp) * 8, destination += bit_size)
	  __pack(reinterpret_cast<const uvalue_type*>(source), reinterpret_cast<uvalue_type*>(destination), bit_size);
    }
    
    template <typename Tp>
    inline
    void unpack(const Tp* source, Tp*  destination, size_t size, size_t bit_size)
    {
      typedef typename __unpack_value<sizeof(Tp)>::value_type uvalue_type;
      
      if (bit_size == sizeof(Tp) * 8) 
	std::copy(source, source + size, destination);
      else
	for (size_t i = 0; i < size; i += sizeof(Tp) * 8, source += bit_size, destination += sizeof(Tp) * 8)
	  __unpack(reinterpret_cast<const uvalue_type*>(source), reinterpret_cast<uvalue_type*>(destination), bit_size);
    }
    
    template <typename Tp>
    inline
    Tp unpack(const Tp* source, size_t pos, size_t bit_size)
    {
      static const size_t shift = utils::bithack::static_floor_log2<sizeof(Tp) * 8>::result;

      if (bit_size == sizeof(Tp) * 8)
	return source[pos];
      else
	return __unpack(source + ((pos >> shift) * bit_size), pos & ((1 << shift) - 1), bit_size);
    }
    
  };
};

#endif
