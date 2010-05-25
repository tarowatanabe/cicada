
#include <stdint.h>

#include <algorithm>
#include <stdexcept>

#include "bitpack.hpp"
#include "bitpack_impl.hpp"
#include "bitpack64_impl.hpp"

namespace utils
{
  namespace bitpack
  {
    
    typedef void (*packer64_type)(const uint64_t* source, uint64_t* destination);
    
    static packer64_type __packer64[63] = {
      __struct_bitpack_impl<uint64_t,8,1>::pack,
      __struct_bitpack_impl<uint64_t,8,2>::pack,
      __struct_bitpack_impl<uint64_t,8,3>::pack,
      __struct_bitpack_impl<uint64_t,8,4>::pack,
      __struct_bitpack_impl<uint64_t,8,5>::pack,
      __struct_bitpack_impl<uint64_t,8,6>::pack,
      __struct_bitpack_impl<uint64_t,8,7>::pack,
      __struct_bitpack_impl<uint64_t,8,8>::pack,
      __struct_bitpack_impl<uint64_t,8,9>::pack,
      __struct_bitpack_impl<uint64_t,8,10>::pack,
      __struct_bitpack_impl<uint64_t,8,11>::pack,
      __struct_bitpack_impl<uint64_t,8,12>::pack,
      __struct_bitpack_impl<uint64_t,8,13>::pack,
      __struct_bitpack_impl<uint64_t,8,14>::pack,
      __struct_bitpack_impl<uint64_t,8,15>::pack,
      __struct_bitpack_impl<uint64_t,8,16>::pack,
      __struct_bitpack_impl<uint64_t,8,17>::pack,
      __struct_bitpack_impl<uint64_t,8,18>::pack,
      __struct_bitpack_impl<uint64_t,8,19>::pack,
      __struct_bitpack_impl<uint64_t,8,20>::pack,
      __struct_bitpack_impl<uint64_t,8,21>::pack,
      __struct_bitpack_impl<uint64_t,8,22>::pack,
      __struct_bitpack_impl<uint64_t,8,23>::pack,
      __struct_bitpack_impl<uint64_t,8,24>::pack,
      __struct_bitpack_impl<uint64_t,8,25>::pack,
      __struct_bitpack_impl<uint64_t,8,26>::pack,
      __struct_bitpack_impl<uint64_t,8,27>::pack,
      __struct_bitpack_impl<uint64_t,8,28>::pack,
      __struct_bitpack_impl<uint64_t,8,29>::pack,
      __struct_bitpack_impl<uint64_t,8,30>::pack,
      __struct_bitpack_impl<uint64_t,8,31>::pack,
      __struct_bitpack_impl<uint64_t,8,32>::pack,
      __struct_bitpack_impl<uint64_t,8,33>::pack,
      __struct_bitpack_impl<uint64_t,8,34>::pack,
      __struct_bitpack_impl<uint64_t,8,35>::pack,
      __struct_bitpack_impl<uint64_t,8,36>::pack,
      __struct_bitpack_impl<uint64_t,8,37>::pack,
      __struct_bitpack_impl<uint64_t,8,38>::pack,
      __struct_bitpack_impl<uint64_t,8,39>::pack,
      __struct_bitpack_impl<uint64_t,8,40>::pack,
      __struct_bitpack_impl<uint64_t,8,41>::pack,
      __struct_bitpack_impl<uint64_t,8,42>::pack,
      __struct_bitpack_impl<uint64_t,8,43>::pack,
      __struct_bitpack_impl<uint64_t,8,44>::pack,
      __struct_bitpack_impl<uint64_t,8,45>::pack,
      __struct_bitpack_impl<uint64_t,8,46>::pack,
      __struct_bitpack_impl<uint64_t,8,47>::pack,
      __struct_bitpack_impl<uint64_t,8,48>::pack,
      __struct_bitpack_impl<uint64_t,8,49>::pack,
      __struct_bitpack_impl<uint64_t,8,50>::pack,
      __struct_bitpack_impl<uint64_t,8,51>::pack,
      __struct_bitpack_impl<uint64_t,8,52>::pack,
      __struct_bitpack_impl<uint64_t,8,53>::pack,
      __struct_bitpack_impl<uint64_t,8,54>::pack,
      __struct_bitpack_impl<uint64_t,8,55>::pack,
      __struct_bitpack_impl<uint64_t,8,56>::pack,
      __struct_bitpack_impl<uint64_t,8,57>::pack,
      __struct_bitpack_impl<uint64_t,8,58>::pack,
      __struct_bitpack_impl<uint64_t,8,59>::pack,
      __struct_bitpack_impl<uint64_t,8,60>::pack,
      __struct_bitpack_impl<uint64_t,8,61>::pack,
      __struct_bitpack_impl<uint64_t,8,62>::pack,
      __struct_bitpack_impl<uint64_t,8,63>::pack,
    };
    static packer64_type __unpacker64[63] = {
      __struct_bitpack_impl<uint64_t,8,1>::unpack,
      __struct_bitpack_impl<uint64_t,8,2>::unpack,
      __struct_bitpack_impl<uint64_t,8,3>::unpack,
      __struct_bitpack_impl<uint64_t,8,4>::unpack,
      __struct_bitpack_impl<uint64_t,8,5>::unpack,
      __struct_bitpack_impl<uint64_t,8,6>::unpack,
      __struct_bitpack_impl<uint64_t,8,7>::unpack,
      __struct_bitpack_impl<uint64_t,8,8>::unpack,
      __struct_bitpack_impl<uint64_t,8,9>::unpack,
      __struct_bitpack_impl<uint64_t,8,10>::unpack,
      __struct_bitpack_impl<uint64_t,8,11>::unpack,
      __struct_bitpack_impl<uint64_t,8,12>::unpack,
      __struct_bitpack_impl<uint64_t,8,13>::unpack,
      __struct_bitpack_impl<uint64_t,8,14>::unpack,
      __struct_bitpack_impl<uint64_t,8,15>::unpack,
      __struct_bitpack_impl<uint64_t,8,16>::unpack,
      __struct_bitpack_impl<uint64_t,8,17>::unpack,
      __struct_bitpack_impl<uint64_t,8,18>::unpack,
      __struct_bitpack_impl<uint64_t,8,19>::unpack,
      __struct_bitpack_impl<uint64_t,8,20>::unpack,
      __struct_bitpack_impl<uint64_t,8,21>::unpack,
      __struct_bitpack_impl<uint64_t,8,22>::unpack,
      __struct_bitpack_impl<uint64_t,8,23>::unpack,
      __struct_bitpack_impl<uint64_t,8,24>::unpack,
      __struct_bitpack_impl<uint64_t,8,25>::unpack,
      __struct_bitpack_impl<uint64_t,8,26>::unpack,
      __struct_bitpack_impl<uint64_t,8,27>::unpack,
      __struct_bitpack_impl<uint64_t,8,28>::unpack,
      __struct_bitpack_impl<uint64_t,8,29>::unpack,
      __struct_bitpack_impl<uint64_t,8,30>::unpack,
      __struct_bitpack_impl<uint64_t,8,31>::unpack,
      __struct_bitpack_impl<uint64_t,8,32>::unpack,
      __struct_bitpack_impl<uint64_t,8,33>::unpack,
      __struct_bitpack_impl<uint64_t,8,34>::unpack,
      __struct_bitpack_impl<uint64_t,8,35>::unpack,
      __struct_bitpack_impl<uint64_t,8,36>::unpack,
      __struct_bitpack_impl<uint64_t,8,37>::unpack,
      __struct_bitpack_impl<uint64_t,8,38>::unpack,
      __struct_bitpack_impl<uint64_t,8,39>::unpack,
      __struct_bitpack_impl<uint64_t,8,40>::unpack,
      __struct_bitpack_impl<uint64_t,8,41>::unpack,
      __struct_bitpack_impl<uint64_t,8,42>::unpack,
      __struct_bitpack_impl<uint64_t,8,43>::unpack,
      __struct_bitpack_impl<uint64_t,8,44>::unpack,
      __struct_bitpack_impl<uint64_t,8,45>::unpack,
      __struct_bitpack_impl<uint64_t,8,46>::unpack,
      __struct_bitpack_impl<uint64_t,8,47>::unpack,
      __struct_bitpack_impl<uint64_t,8,48>::unpack,
      __struct_bitpack_impl<uint64_t,8,49>::unpack,
      __struct_bitpack_impl<uint64_t,8,50>::unpack,
      __struct_bitpack_impl<uint64_t,8,51>::unpack,
      __struct_bitpack_impl<uint64_t,8,52>::unpack,
      __struct_bitpack_impl<uint64_t,8,53>::unpack,
      __struct_bitpack_impl<uint64_t,8,54>::unpack,
      __struct_bitpack_impl<uint64_t,8,55>::unpack,
      __struct_bitpack_impl<uint64_t,8,56>::unpack,
      __struct_bitpack_impl<uint64_t,8,57>::unpack,
      __struct_bitpack_impl<uint64_t,8,58>::unpack,
      __struct_bitpack_impl<uint64_t,8,59>::unpack,
      __struct_bitpack_impl<uint64_t,8,60>::unpack,
      __struct_bitpack_impl<uint64_t,8,61>::unpack,
      __struct_bitpack_impl<uint64_t,8,62>::unpack,
      __struct_bitpack_impl<uint64_t,8,63>::unpack,
    };
    
    void __pack(const uint64_t* source, uint64_t* destination, size_t bit_size)
    {
      __packer64[bit_size - 1](source, destination);
    }

    void __unpack(const uint64_t* source, uint64_t* destination, size_t bit_size)
    {
      __unpacker64[bit_size - 1](source, destination);
    }    
  };  
};
