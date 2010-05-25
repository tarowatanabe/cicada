
#include <iostream>
#include <vector>

#include "pfor_vector.hpp"
#include "pfor_device.hpp"

template <typename Tp>
struct Task
{

  void operator()()
  {
    std::ostringstream stream_value_size;
    stream_value_size << (sizeof(Tp) * 8);
    
    const std::string value_size = stream_value_size.str();

    for (int i = 0; i < 16; ++ i) {
      std::cerr << "iteration: " << i << std::endl;
    
      std::vector<Tp> integers;
      for (int j = 0; j < 256 * 256 + 10; ++ j)
	integers.push_back(random() % (1024 * 64));
      
      utils::pfor_vector<Tp> pfor_ints;
      boost::iostreams::filtering_ostream os;
      os.push(utils::pfor_sink<Tp>(std::string("tmptmp.pfor.iostream.") + value_size));
      
      for (int j = 0; j < 256 * 256 + 10; ++ j) {
	pfor_ints.push_back(integers[j]);
	os.write((char*) &(integers[j]), sizeof(Tp));
      }
      pfor_ints.build();
      pfor_ints.write(std::string("tmptmp.pfor.") + value_size);
      os.pop();
      
      utils::pfor_vector_mapped<Tp> pfor_ints_mapped(std::string("tmptmp.pfor.") + value_size);
      utils::pfor_vector_mapped<Tp> pfor_ints_stream(std::string("tmptmp.pfor.iostream.") + value_size);
      for (int j = 0; j < 256 * 256 + 10; ++ j) {
	
	const Tp value = pfor_ints[j];
	const Tp value_mapped = pfor_ints_mapped[j];
	const Tp value_stream = pfor_ints_stream[j];
	
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

int main(int argc, char** argv)
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
}
