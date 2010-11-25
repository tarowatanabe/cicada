//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>

#include "packed_vector.hpp"
#include "packed_device.hpp"

template <typename Tp>
struct Task
{

  void operator()()
  {
    for (int i = 0; i < 16; ++ i) {
      std::cerr << "iteration: " << i << std::endl;
    
      std::vector<Tp> integers;
      for (int j = 0; j < 256 * 256 + 10; ++ j)
	integers.push_back(random() % (1024 * 64));
      
      utils::packed_vector<Tp> packed_ints;
      boost::iostreams::filtering_ostream os;
      os.push(utils::packed_sink<Tp>("tmptmp.packed.iostream"));
      
      for (int j = 0; j < 256 * 256 + 10; ++ j) {
	packed_ints.push_back(integers[j]);
	os.write((char*) &(integers[j]), sizeof(Tp));
      }
      packed_ints.build();
      packed_ints.write("tmptmp.packed");
      os.pop();
      
      utils::packed_vector_mapped<Tp> packed_ints_mapped("tmptmp.packed");
      utils::packed_vector_mapped<Tp> packed_ints_stream("tmptmp.packed.iostream");
      for (int j = 0; j < 256 * 256 + 10; ++ j) {
	
	const Tp value = packed_ints[j];
	const Tp value_mapped = packed_ints_mapped[j];
	const Tp value_stream = packed_ints_stream[j];
	
	if (value != integers[j])
	  std::cerr << "DIFFER(raw   ): i = " << j << " " << int64_t(value) << " " << int64_t(integers[j]) << std::endl;
	
	if (value_mapped != integers[j])
	  std::cerr << "DIFFER(mapped): i = " << j << " " << int64_t(value_mapped) << " " << int64_t(integers[j]) << std::endl;
	
	if (value_stream != integers[j])
	  std::cerr << "DIFFER(stream): i = " << j << " " << int64_t(value_stream) << " " << int64_t(integers[j]) << std::endl;
      }
    }
  }
};

int main(int argc, char**argv)
{
  srandom(time(0) * getpid());
  
  std::cout << "64bit" << std::endl;
  Task<int64_t> task_int64_t;
  task_int64_t();

  std::cout << "32bit" << std::endl;
  Task<int> task_int;
  task_int();
  
  std::cout << "16bit" << std::endl;
  Task<int16_t> task_int16_t;
  task_int16_t();
  
  std::cout << "8bit" << std::endl;
  Task<char> task_char;
  task_char();
}
