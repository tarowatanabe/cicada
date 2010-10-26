
#include <iostream>

#include <climits>
#include <vector>

#include "group_aligned_code.hpp"

#include <boost/numeric/conversion/bounds.hpp>


int main(int argc, char** argv)
{
  char buffer[32];

  uint32_t max_int32 = UINT_MAX;
  uint32_t max_int32_decoded;
  utils::group_aligned_encode(max_int32, buffer, 0);
  utils::group_aligned_decode(max_int32_decoded, buffer, 0);
    
  if (max_int32 != max_int32_decoded)
    std::cerr << "MAX UINT32 differ!" << std::endl;

  srandom(time(0) * getpid());
  
  for (int i = 0; i < 1024 * 64; ++ i) {
    const uint32_t value = random();
    uint32_t value_decoded;

    const uint32_t pos = random() & 0x03;
    
    const size_t encode_offset = utils::group_aligned_encode(value, buffer, pos);
    const size_t decode_offset = utils::group_aligned_decode(value_decoded, buffer, pos);
    
    if (value != value_decoded)
      std::cerr << "DIFFER for 32-bit" << std::endl;
    
    if (encode_offset != decode_offset)
      std::cerr << "DIFFER for 62-bit" << std::endl;
  }

  for (int i = 0; i < 1024; ++ i) {
    const size_t num_sample = random() & (1024 - 1);

    std::vector<char> buffer(num_sample * 8);
    std::vector<int>  values;
    {
      std::vector<char>::iterator biter = buffer.begin();
      std::vector<char>::iterator hiter = buffer.begin();
      for (size_t j = 0; j < num_sample; ++ j) {
	const uint32_t value = random();
	
	const size_t offset = utils::group_aligned_encode(value, &(*hiter), j);
	
	biter = hiter + offset;
	hiter += offset & (- size_t((j & 0x03) == 0x03));
	
	values.push_back(value);
      }
      buffer.resize(biter - buffer.begin());
    }
    
    std::vector<char>::iterator biter = buffer.begin();
    std::vector<char>::iterator hiter = buffer.begin();
    
    size_t decoded_size = 0;
    for (int j = 0; biter != buffer.end(); ++ j, ++ decoded_size) {
      
      uint32_t value;
      const size_t offset = utils::group_aligned_decode(value, &(*hiter), j);
      
      biter = hiter + offset;
      hiter += offset & (- size_t((j & 0x03) == 0x03));

      if (value != values[j])
	std::cerr << "differ..." << std::endl;
    }

    if (decoded_size != num_sample)
      std::cerr << "different sample size..." << std::endl;
  }

}
