//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdint.h>

#include <algorithm>
#include <stdexcept>

#include "bitpack.hpp"
#include "bitpack_impl.hpp"
#include "bitpack32_impl.hpp"

namespace utils
{
  namespace bitpack
  {
    
    typedef void (*packer32_type)(const uint32_t* source, uint32_t* destination);

    static packer32_type __packer32[31] = {
      __struct_bitpack_impl<uint32_t,4,1>::pack,
      __struct_bitpack_impl<uint32_t,4,2>::pack,
      __struct_bitpack_impl<uint32_t,4,3>::pack,
      __struct_bitpack_impl<uint32_t,4,4>::pack,
      __struct_bitpack_impl<uint32_t,4,5>::pack,
      __struct_bitpack_impl<uint32_t,4,6>::pack,
      __struct_bitpack_impl<uint32_t,4,7>::pack,
      __struct_bitpack_impl<uint32_t,4,8>::pack,
      __struct_bitpack_impl<uint32_t,4,9>::pack,
      __struct_bitpack_impl<uint32_t,4,10>::pack,
      __struct_bitpack_impl<uint32_t,4,11>::pack,
      __struct_bitpack_impl<uint32_t,4,12>::pack,
      __struct_bitpack_impl<uint32_t,4,13>::pack,
      __struct_bitpack_impl<uint32_t,4,14>::pack,
      __struct_bitpack_impl<uint32_t,4,15>::pack,
      __struct_bitpack_impl<uint32_t,4,16>::pack,
      __struct_bitpack_impl<uint32_t,4,17>::pack,
      __struct_bitpack_impl<uint32_t,4,18>::pack,
      __struct_bitpack_impl<uint32_t,4,19>::pack,
      __struct_bitpack_impl<uint32_t,4,20>::pack,
      __struct_bitpack_impl<uint32_t,4,21>::pack,
      __struct_bitpack_impl<uint32_t,4,22>::pack,
      __struct_bitpack_impl<uint32_t,4,23>::pack,
      __struct_bitpack_impl<uint32_t,4,24>::pack,
      __struct_bitpack_impl<uint32_t,4,25>::pack,
      __struct_bitpack_impl<uint32_t,4,26>::pack,
      __struct_bitpack_impl<uint32_t,4,27>::pack,
      __struct_bitpack_impl<uint32_t,4,28>::pack,
      __struct_bitpack_impl<uint32_t,4,29>::pack,
      __struct_bitpack_impl<uint32_t,4,30>::pack,
      __struct_bitpack_impl<uint32_t,4,31>::pack,
    };
    static packer32_type __unpacker32[31] = {
      __struct_bitpack_impl<uint32_t,4,1>::unpack,
      __struct_bitpack_impl<uint32_t,4,2>::unpack,
      __struct_bitpack_impl<uint32_t,4,3>::unpack,
      __struct_bitpack_impl<uint32_t,4,4>::unpack,
      __struct_bitpack_impl<uint32_t,4,5>::unpack,
      __struct_bitpack_impl<uint32_t,4,6>::unpack,
      __struct_bitpack_impl<uint32_t,4,7>::unpack,
      __struct_bitpack_impl<uint32_t,4,8>::unpack,
      __struct_bitpack_impl<uint32_t,4,9>::unpack,
      __struct_bitpack_impl<uint32_t,4,10>::unpack,
      __struct_bitpack_impl<uint32_t,4,11>::unpack,
      __struct_bitpack_impl<uint32_t,4,12>::unpack,
      __struct_bitpack_impl<uint32_t,4,13>::unpack,
      __struct_bitpack_impl<uint32_t,4,14>::unpack,
      __struct_bitpack_impl<uint32_t,4,15>::unpack,
      __struct_bitpack_impl<uint32_t,4,16>::unpack,
      __struct_bitpack_impl<uint32_t,4,17>::unpack,
      __struct_bitpack_impl<uint32_t,4,18>::unpack,
      __struct_bitpack_impl<uint32_t,4,19>::unpack,
      __struct_bitpack_impl<uint32_t,4,20>::unpack,
      __struct_bitpack_impl<uint32_t,4,21>::unpack,
      __struct_bitpack_impl<uint32_t,4,22>::unpack,
      __struct_bitpack_impl<uint32_t,4,23>::unpack,
      __struct_bitpack_impl<uint32_t,4,24>::unpack,
      __struct_bitpack_impl<uint32_t,4,25>::unpack,
      __struct_bitpack_impl<uint32_t,4,26>::unpack,
      __struct_bitpack_impl<uint32_t,4,27>::unpack,
      __struct_bitpack_impl<uint32_t,4,28>::unpack,
      __struct_bitpack_impl<uint32_t,4,29>::unpack,
      __struct_bitpack_impl<uint32_t,4,30>::unpack,
      __struct_bitpack_impl<uint32_t,4,31>::unpack,
    };    

    void __pack(const uint32_t* source, uint32_t* destination, size_t bit_size)
    {
      __packer32[bit_size - 1](source, destination);
    }
    
    void __unpack(const uint32_t* source, uint32_t* destination, size_t bit_size)
    {
      __unpacker32[bit_size - 1](source, destination);
    }
  };  
};
