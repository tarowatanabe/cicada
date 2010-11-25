//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <climits>

#include "byte_aligned_code.hpp"

#include <boost/numeric/conversion/bounds.hpp>


int main(int argc, char** argv)
{
  char buffer[16];

  uint32_t max_int32 = UINT_MAX;
  uint32_t max_int32_decoded;
  utils::byte_aligned_encode(max_int32, buffer);
  utils::byte_aligned_decode(max_int32_decoded, buffer);
    
  if (max_int32 != max_int32_decoded)
    std::cerr << "MAX UINT32 differ!" << std::endl;
  
  srandom(time(0) * getpid());
  
  for (int i = 0; i < 1024; ++ i) {
    const uint32_t value = random();
    uint32_t value_decoded;

    const size_t encode_size = utils::byte_aligned_encode(value, buffer);
    const size_t decode_size = utils::byte_aligned_decode(value_decoded, buffer);
    
    if (value != value_decoded)
      std::cerr << "DIFFER for 32-bit" << std::endl;
    
    if (encode_size != decode_size)
      std::cerr << "DIFFER for 62-bit" << std::endl;
  }

  uint64_t max_int64 = boost::numeric::bounds<uint64_t>::highest();
  uint64_t max_int64_decoded;
  utils::byte_aligned_encode(max_int64, buffer);
  utils::byte_aligned_decode(max_int64_decoded, buffer);
    
  if (max_int64 != max_int64_decoded)
    std::cerr << "MAX UINT64 differ!" << std::endl;

  for (int i = 0; i < 1024; ++ i) {
    const uint64_t value = (uint64_t(random()) << 32) | random();
    uint64_t value_decoded;

    const size_t encode_size = utils::byte_aligned_encode(value, buffer);
    const size_t decode_size = utils::byte_aligned_decode(value_decoded, buffer);
    
    if (value != value_decoded)
      std::cerr << "DIFFER for 32-bit" << std::endl;
    
    if (encode_size != decode_size)
      std::cerr << "DIFFER for 64-bit" << std::endl;
  }
      
}
