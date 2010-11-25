//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdint.h>

#include <algorithm>
#include <stdexcept>

#include "bitpack.hpp"
#include "bitpack_impl.hpp"
#include "bitpack16_impl.hpp"

namespace utils
{
  namespace bitpack
  {
    
    typedef void (*packer16_type)(const uint16_t* source, uint16_t* destination);
    
    static packer16_type __packer16[15] = {
      __struct_bitpack_impl<uint16_t,2,1>::pack,
      __struct_bitpack_impl<uint16_t,2,2>::pack,
      __struct_bitpack_impl<uint16_t,2,3>::pack,
      __struct_bitpack_impl<uint16_t,2,4>::pack,
      __struct_bitpack_impl<uint16_t,2,5>::pack,
      __struct_bitpack_impl<uint16_t,2,6>::pack,
      __struct_bitpack_impl<uint16_t,2,7>::pack,
      __struct_bitpack_impl<uint16_t,2,8>::pack,
      __struct_bitpack_impl<uint16_t,2,9>::pack,
      __struct_bitpack_impl<uint16_t,2,10>::pack,
      __struct_bitpack_impl<uint16_t,2,11>::pack,
      __struct_bitpack_impl<uint16_t,2,12>::pack,
      __struct_bitpack_impl<uint16_t,2,13>::pack,
      __struct_bitpack_impl<uint16_t,2,14>::pack,
      __struct_bitpack_impl<uint16_t,2,15>::pack,
    };
    static packer16_type __unpacker16[15] = {
      __struct_bitpack_impl<uint16_t,2,1>::unpack,
      __struct_bitpack_impl<uint16_t,2,2>::unpack,
      __struct_bitpack_impl<uint16_t,2,3>::unpack,
      __struct_bitpack_impl<uint16_t,2,4>::unpack,
      __struct_bitpack_impl<uint16_t,2,5>::unpack,
      __struct_bitpack_impl<uint16_t,2,6>::unpack,
      __struct_bitpack_impl<uint16_t,2,7>::unpack,
      __struct_bitpack_impl<uint16_t,2,8>::unpack,
      __struct_bitpack_impl<uint16_t,2,9>::unpack,
      __struct_bitpack_impl<uint16_t,2,10>::unpack,
      __struct_bitpack_impl<uint16_t,2,11>::unpack,
      __struct_bitpack_impl<uint16_t,2,12>::unpack,
      __struct_bitpack_impl<uint16_t,2,13>::unpack,
      __struct_bitpack_impl<uint16_t,2,14>::unpack,
      __struct_bitpack_impl<uint16_t,2,15>::unpack,
    };
    
    void __pack(const uint16_t* source, uint16_t* destination, size_t bit_size)
    {
      __packer16[bit_size - 1](source, destination);
    }
    
    void __unpack(const uint16_t* source, uint16_t* destination, size_t bit_size)
    {
      __unpacker16[bit_size - 1](source, destination);
    }
  };  
};
