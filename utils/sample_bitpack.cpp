//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>

#include <utils/bitpack.hpp>

template <typename Tp>
struct Task
{
  
  void operator()()
  {
    std::vector<Tp> vec;
    std::vector<Tp> packed;
    std::vector<Tp> decoded;

    const size_t value_size = sizeof(Tp) * 8;
    
    for (size_t bits = 1; bits <= value_size; ++ bits) {    
      vec.clear();
      vec.resize(value_size * 2, 0);
      decoded.resize(value_size * 2);      
      packed.resize(bits * 2);
      
      for (size_t i = 0; i < vec.size(); ++ i)
	vec[i] = random() & ((1 << bits) - 1);
      
      std::cerr << "# of bits: " << bits << std::endl;
      
      std::cerr << "Encoding" << std::endl;
      utils::bitpack::pack(&(*vec.begin()), &(*packed.begin()), value_size * 2, bits);
      
      std::cerr << "Decoding" << std::endl;
      utils::bitpack::unpack(&(*packed.begin()), &(*decoded.begin()), value_size * 2, bits);
    
    
      for (size_t i = 0; i < vec.size(); ++ i) {
	if (vec[i] != decoded[i])
	  std::cerr << "i = " << i << " orig: " << vec[i] << " decoded: " << decoded[i] << std::endl;
	
	if (vec[i] != utils::bitpack::unpack(&(*packed.begin()), i, bits))
	  std::cerr << "i = " << i << " orig: " << vec[i] << " access: " << utils::bitpack::unpack(&(*packed.begin()), i, bits) << std::endl;
      }
    }
  }
};


int main(int argc, char** argv)
{
  

  srandom(time(0) * getpid());
  
  Task<int> task_int;
  task_int();
  
  Task<int64_t> task_64;
  task_64();
  
  Task<int16_t> task_16;
  task_16();
}
