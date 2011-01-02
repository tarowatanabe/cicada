//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdint.h>

#include <algorithm>
#include <stdexcept>

#include "bitpack.hpp"
#include "bitpack_impl.hpp"
#include "bitpack8_impl.hpp"

namespace utils
{
  namespace bitpack
  {
    typedef void (*packer8_type)(const uint8_t* source, uint8_t* destination);
    
    static packer8_type  __packer8[7] = {
      __struct_bitpack_impl<uint8_t,1,1>::pack,
      __struct_bitpack_impl<uint8_t,1,2>::pack,
      __struct_bitpack_impl<uint8_t,1,3>::pack,
      __struct_bitpack_impl<uint8_t,1,4>::pack,
      __struct_bitpack_impl<uint8_t,1,5>::pack,
      __struct_bitpack_impl<uint8_t,1,6>::pack,
      __struct_bitpack_impl<uint8_t,1,7>::pack,
    };
    static packer8_type  __unpacker8[7] = {
      __struct_bitpack_impl<uint8_t,1,1>::unpack,
      __struct_bitpack_impl<uint8_t,1,2>::unpack,
      __struct_bitpack_impl<uint8_t,1,3>::unpack,
      __struct_bitpack_impl<uint8_t,1,4>::unpack,
      __struct_bitpack_impl<uint8_t,1,5>::unpack,
      __struct_bitpack_impl<uint8_t,1,6>::unpack,
      __struct_bitpack_impl<uint8_t,1,7>::unpack,
    };
    
    void __pack(const uint8_t*  source, uint8_t*  destination, size_t bit_size)
    {
      __packer8[bit_size - 1](source, destination);
    }
    
    void __unpack(const uint8_t*  source, uint8_t*  destination, size_t bit_size)
    {
      __unpacker8[bit_size - 1](source, destination);
    }
    
  };  
};
