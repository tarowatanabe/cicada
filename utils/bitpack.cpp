//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdint.h>

#include <algorithm>
#include <stdexcept>

#include "bitpack.hpp"

namespace utils
{
  namespace bitpack
  {
    
    typedef void (*packer8_type)(const uint8_t* source, uint8_t* destination);
    typedef void (*packer16_type)(const uint16_t* source, uint16_t* destination);
    typedef void (*packer32_type)(const uint32_t* source, uint32_t* destination);
    typedef void (*packer64_type)(const uint64_t* source, uint64_t* destination);
    
    static packer8_type  __packer8[7] = {
      __struct_bitpack_loop<uint8_t,1,1,8>::pack,
      __struct_bitpack_loop<uint8_t,1,2,8>::pack,
      __struct_bitpack_loop<uint8_t,1,3,8>::pack,
      __struct_bitpack_loop<uint8_t,1,4,8>::pack,
      __struct_bitpack_loop<uint8_t,1,5,8>::pack,
      __struct_bitpack_loop<uint8_t,1,6,8>::pack,
      __struct_bitpack_loop<uint8_t,1,7,8>::pack,
    };
    static packer8_type  __unpacker8[7] = {
      __struct_bitpack_loop<uint8_t,1,1,8>::unpack,
      __struct_bitpack_loop<uint8_t,1,2,8>::unpack,
      __struct_bitpack_loop<uint8_t,1,3,8>::unpack,
      __struct_bitpack_loop<uint8_t,1,4,8>::unpack,
      __struct_bitpack_loop<uint8_t,1,5,8>::unpack,
      __struct_bitpack_loop<uint8_t,1,6,8>::unpack,
      __struct_bitpack_loop<uint8_t,1,7,8>::unpack,
    };
    static packer16_type __packer16[15] = {
      __struct_bitpack_loop<uint16_t,2,1,16>::pack,
      __struct_bitpack_loop<uint16_t,2,2,16>::pack,
      __struct_bitpack_loop<uint16_t,2,3,16>::pack,
      __struct_bitpack_loop<uint16_t,2,4,16>::pack,
      __struct_bitpack_loop<uint16_t,2,5,16>::pack,
      __struct_bitpack_loop<uint16_t,2,6,16>::pack,
      __struct_bitpack_loop<uint16_t,2,7,16>::pack,
      __struct_bitpack_loop<uint16_t,2,8,16>::pack,
      __struct_bitpack_loop<uint16_t,2,9,16>::pack,
      __struct_bitpack_loop<uint16_t,2,10,16>::pack,
      __struct_bitpack_loop<uint16_t,2,11,16>::pack,
      __struct_bitpack_loop<uint16_t,2,12,16>::pack,
      __struct_bitpack_loop<uint16_t,2,13,16>::pack,
      __struct_bitpack_loop<uint16_t,2,14,16>::pack,
      __struct_bitpack_loop<uint16_t,2,15,16>::pack,
    };
    static packer16_type __unpacker16[15] = {
      __struct_bitpack_loop<uint16_t,2,1,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,2,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,3,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,4,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,5,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,6,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,7,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,8,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,9,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,10,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,11,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,12,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,13,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,14,16>::unpack,
      __struct_bitpack_loop<uint16_t,2,15,16>::unpack,
    };
    static packer32_type __packer32[31] = {
      __struct_bitpack_loop<uint32_t,4,1,32>::pack,
      __struct_bitpack_loop<uint32_t,4,2,32>::pack,
      __struct_bitpack_loop<uint32_t,4,3,32>::pack,
      __struct_bitpack_loop<uint32_t,4,4,32>::pack,
      __struct_bitpack_loop<uint32_t,4,5,32>::pack,
      __struct_bitpack_loop<uint32_t,4,6,32>::pack,
      __struct_bitpack_loop<uint32_t,4,7,32>::pack,
      __struct_bitpack_loop<uint32_t,4,8,32>::pack,
      __struct_bitpack_loop<uint32_t,4,9,32>::pack,
      __struct_bitpack_loop<uint32_t,4,10,32>::pack,
      __struct_bitpack_loop<uint32_t,4,11,32>::pack,
      __struct_bitpack_loop<uint32_t,4,12,32>::pack,
      __struct_bitpack_loop<uint32_t,4,13,32>::pack,
      __struct_bitpack_loop<uint32_t,4,14,32>::pack,
      __struct_bitpack_loop<uint32_t,4,15,32>::pack,
      __struct_bitpack_loop<uint32_t,4,16,32>::pack,
      __struct_bitpack_loop<uint32_t,4,17,32>::pack,
      __struct_bitpack_loop<uint32_t,4,18,32>::pack,
      __struct_bitpack_loop<uint32_t,4,19,32>::pack,
      __struct_bitpack_loop<uint32_t,4,20,32>::pack,
      __struct_bitpack_loop<uint32_t,4,21,32>::pack,
      __struct_bitpack_loop<uint32_t,4,22,32>::pack,
      __struct_bitpack_loop<uint32_t,4,23,32>::pack,
      __struct_bitpack_loop<uint32_t,4,24,32>::pack,
      __struct_bitpack_loop<uint32_t,4,25,32>::pack,
      __struct_bitpack_loop<uint32_t,4,26,32>::pack,
      __struct_bitpack_loop<uint32_t,4,27,32>::pack,
      __struct_bitpack_loop<uint32_t,4,28,32>::pack,
      __struct_bitpack_loop<uint32_t,4,29,32>::pack,
      __struct_bitpack_loop<uint32_t,4,30,32>::pack,
      __struct_bitpack_loop<uint32_t,4,31,32>::pack,
    };
    static packer32_type __unpacker32[31] = {
      __struct_bitpack_loop<uint32_t,4,1,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,2,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,3,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,4,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,5,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,6,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,7,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,8,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,9,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,10,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,11,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,12,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,13,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,14,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,15,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,16,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,17,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,18,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,19,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,20,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,21,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,22,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,23,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,24,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,25,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,26,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,27,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,28,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,29,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,30,32>::unpack,
      __struct_bitpack_loop<uint32_t,4,31,32>::unpack,
    };
    static packer64_type __packer64[63] = {
      __struct_bitpack_loop<uint64_t,8,1,64>::pack,
      __struct_bitpack_loop<uint64_t,8,2,64>::pack,
      __struct_bitpack_loop<uint64_t,8,3,64>::pack,
      __struct_bitpack_loop<uint64_t,8,4,64>::pack,
      __struct_bitpack_loop<uint64_t,8,5,64>::pack,
      __struct_bitpack_loop<uint64_t,8,6,64>::pack,
      __struct_bitpack_loop<uint64_t,8,7,64>::pack,
      __struct_bitpack_loop<uint64_t,8,8,64>::pack,
      __struct_bitpack_loop<uint64_t,8,9,64>::pack,
      __struct_bitpack_loop<uint64_t,8,10,64>::pack,
      __struct_bitpack_loop<uint64_t,8,11,64>::pack,
      __struct_bitpack_loop<uint64_t,8,12,64>::pack,
      __struct_bitpack_loop<uint64_t,8,13,64>::pack,
      __struct_bitpack_loop<uint64_t,8,14,64>::pack,
      __struct_bitpack_loop<uint64_t,8,15,64>::pack,
      __struct_bitpack_loop<uint64_t,8,16,64>::pack,
      __struct_bitpack_loop<uint64_t,8,17,64>::pack,
      __struct_bitpack_loop<uint64_t,8,18,64>::pack,
      __struct_bitpack_loop<uint64_t,8,19,64>::pack,
      __struct_bitpack_loop<uint64_t,8,20,64>::pack,
      __struct_bitpack_loop<uint64_t,8,21,64>::pack,
      __struct_bitpack_loop<uint64_t,8,22,64>::pack,
      __struct_bitpack_loop<uint64_t,8,23,64>::pack,
      __struct_bitpack_loop<uint64_t,8,24,64>::pack,
      __struct_bitpack_loop<uint64_t,8,25,64>::pack,
      __struct_bitpack_loop<uint64_t,8,26,64>::pack,
      __struct_bitpack_loop<uint64_t,8,27,64>::pack,
      __struct_bitpack_loop<uint64_t,8,28,64>::pack,
      __struct_bitpack_loop<uint64_t,8,29,64>::pack,
      __struct_bitpack_loop<uint64_t,8,30,64>::pack,
      __struct_bitpack_loop<uint64_t,8,31,64>::pack,
      __struct_bitpack_loop<uint64_t,8,32,64>::pack,
      __struct_bitpack_loop<uint64_t,8,33,64>::pack,
      __struct_bitpack_loop<uint64_t,8,34,64>::pack,
      __struct_bitpack_loop<uint64_t,8,35,64>::pack,
      __struct_bitpack_loop<uint64_t,8,36,64>::pack,
      __struct_bitpack_loop<uint64_t,8,37,64>::pack,
      __struct_bitpack_loop<uint64_t,8,38,64>::pack,
      __struct_bitpack_loop<uint64_t,8,39,64>::pack,
      __struct_bitpack_loop<uint64_t,8,40,64>::pack,
      __struct_bitpack_loop<uint64_t,8,41,64>::pack,
      __struct_bitpack_loop<uint64_t,8,42,64>::pack,
      __struct_bitpack_loop<uint64_t,8,43,64>::pack,
      __struct_bitpack_loop<uint64_t,8,44,64>::pack,
      __struct_bitpack_loop<uint64_t,8,45,64>::pack,
      __struct_bitpack_loop<uint64_t,8,46,64>::pack,
      __struct_bitpack_loop<uint64_t,8,47,64>::pack,
      __struct_bitpack_loop<uint64_t,8,48,64>::pack,
      __struct_bitpack_loop<uint64_t,8,49,64>::pack,
      __struct_bitpack_loop<uint64_t,8,50,64>::pack,
      __struct_bitpack_loop<uint64_t,8,51,64>::pack,
      __struct_bitpack_loop<uint64_t,8,52,64>::pack,
      __struct_bitpack_loop<uint64_t,8,53,64>::pack,
      __struct_bitpack_loop<uint64_t,8,54,64>::pack,
      __struct_bitpack_loop<uint64_t,8,55,64>::pack,
      __struct_bitpack_loop<uint64_t,8,56,64>::pack,
      __struct_bitpack_loop<uint64_t,8,57,64>::pack,
      __struct_bitpack_loop<uint64_t,8,58,64>::pack,
      __struct_bitpack_loop<uint64_t,8,59,64>::pack,
      __struct_bitpack_loop<uint64_t,8,60,64>::pack,
      __struct_bitpack_loop<uint64_t,8,61,64>::pack,
      __struct_bitpack_loop<uint64_t,8,62,64>::pack,
      __struct_bitpack_loop<uint64_t,8,63,64>::pack,
    };
    static packer64_type __unpacker64[63] = {
      __struct_bitpack_loop<uint64_t,8,1,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,2,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,3,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,4,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,5,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,6,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,7,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,8,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,9,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,10,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,11,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,12,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,13,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,14,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,15,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,16,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,17,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,18,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,19,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,20,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,21,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,22,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,23,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,24,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,25,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,26,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,27,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,28,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,29,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,30,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,31,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,32,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,33,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,34,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,35,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,36,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,37,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,38,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,39,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,40,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,41,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,42,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,43,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,44,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,45,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,46,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,47,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,48,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,49,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,50,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,51,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,52,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,53,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,54,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,55,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,56,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,57,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,58,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,59,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,60,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,61,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,62,64>::unpack,
      __struct_bitpack_loop<uint64_t,8,63,64>::unpack,
    };
    
    void __pack(const uint8_t*  source, uint8_t*  destination, size_t bit_size)
    {
      __packer8[bit_size - 1](source, destination);
    }
    void __pack(const uint16_t* source, uint16_t* destination, size_t bit_size)
    {
      __packer16[bit_size - 1](source, destination);
    }
    void __pack(const uint32_t* source, uint32_t* destination, size_t bit_size)
    {
      __packer32[bit_size - 1](source, destination);
    }
    void __pack(const uint64_t* source, uint64_t* destination, size_t bit_size)
    {
      __packer64[bit_size - 1](source, destination);
    }
    
    void __unpack(const uint8_t*  source, uint8_t*  destination, size_t bit_size)
    {
      __unpacker8[bit_size - 1](source, destination);
    }
    void __unpack(const uint16_t* source, uint16_t* destination, size_t bit_size)
    {
      __unpacker16[bit_size - 1](source, destination);
    }
    void __unpack(const uint32_t* source, uint32_t* destination, size_t bit_size)
    {
      __unpacker32[bit_size - 1](source, destination);
    }
    void __unpack(const uint64_t* source, uint64_t* destination, size_t bit_size)
    {
      __unpacker64[bit_size - 1](source, destination);
    }
  };  
};
