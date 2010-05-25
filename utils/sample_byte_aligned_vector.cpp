
#include <unistd.h>
#include <time.h>

#include <iostream>
#include <vector>
#include <algorithm>

#include "utils/byte_aligned_vector.hpp"
#include "utils/byte_aligned_delta_vector.hpp"

#include "utils/byte_aligned_pair_delta_vector.hpp"

template <typename Packed, typename Vector>
struct Task
{
  static inline
  void run()
  {
    typedef Packed packed_vector_type;
    typedef Vector vector_type;
    
    for (int iter = 0; iter < 10; ++ iter) {
      vector_type vec;
    
      for (int i = 0; i < 1024 * 128; ++ i)
	vec.push_back(random() - random());

      std::sort(vec.begin(), vec.end());
    
      packed_vector_type packed(vec.begin(), vec.end());

      std::cerr << "size: " << (vec.size() * sizeof(typename vector_type::value_type)) 
		<< " compressed: " << packed.compressed_size()
		<< std::endl;
    
      size_t size_vec = 0;
      {
	typename packed_vector_type::const_iterator piter = packed.begin();
	typename vector_type::const_iterator iter = vec.begin();
	for (/**/; iter != vec.end(); ++ iter, ++ piter, ++ size_vec)
	  if (*iter != *piter)
	    std::cerr << "differ: " << *iter << " " << *piter << std::endl;
      }
    

      size_t size_packed = 0;
      {
	typename packed_vector_type::const_iterator piter = packed.begin();
	typename vector_type::const_iterator iter = vec.begin();
	for (/**/; piter != packed.end(); ++ iter, ++ piter, ++ size_packed)
	  if (*iter != *piter)
	    std::cerr << "differ: " << *iter << " " << *piter << std::endl;
      }
    
      if (size_vec != size_packed)
	std::cerr << "differ: " << size_vec << " " << size_packed << std::endl;
    
    }
    
  }
  
};

template <typename Packed, typename Vector>
struct TaskPair
{
  static inline
  void run()
  {
    typedef Packed packed_vector_type;
    typedef Vector vector_type;
    
    for (int iter = 0; iter < 10; ++ iter) {
      vector_type vec;
    
      for (int i = 0; i < 1024 * 128; ++ i)
	vec.push_back(std::make_pair(i % 128, (random() % 1024) * (random() & 0x01 ? -1 : 1)));

      std::sort(vec.begin(), vec.end());
    
      packed_vector_type packed(vec.begin(), vec.end());

      std::cerr << "size: " << (vec.size() * sizeof(typename vector_type::value_type)) 
		<< " compressed: " << packed.compressed_size()
		<< std::endl;
    
      size_t size_vec = 0;
      {
	typename packed_vector_type::const_iterator piter = packed.begin();
	typename vector_type::const_iterator iter = vec.begin();
	for (/**/; iter != vec.end(); ++ iter, ++ piter, ++ size_vec)
	  if (*iter != *piter)
	    std::cerr << "differ: " << iter->first << ' ' << iter->second
		      << " " << piter->first << ' ' << piter->second << std::endl;
      }
    

      size_t size_packed = 0;
      {
	typename packed_vector_type::const_iterator piter = packed.begin();
	typename vector_type::const_iterator iter = vec.begin();
	for (/**/; piter != packed.end(); ++ iter, ++ piter, ++ size_packed)
	  if (*iter != *piter)
	    std::cerr << "differ: " << iter->first << ' ' << iter->second
		      << " " << piter->first << ' ' << piter->second << std::endl;
      }
    
      if (size_vec != size_packed)
	std::cerr << "differ: " << size_vec << " " << size_packed << std::endl;
    }
    
  }
  
};

int main(int argc, char** argv)
{
  srandom(time(0) * getpid());
  
  std::cerr << "aligned vector" << std::endl;
  Task<utils::byte_aligned_vector<int>, std::vector<int> >::run();
  
  std::cerr << "aligned delta vector" << std::endl;
  Task<utils::byte_aligned_delta_vector<int>, std::vector<int> >::run();
  
  std::cerr << "aligned pair delta vector" << std::endl;
  TaskPair<utils::byte_aligned_pair_delta_vector<std::pair<int, int> >, std::vector<std::pair<int, int> > >::run();
  
}
