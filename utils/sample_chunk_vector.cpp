//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <string>
#include <cstdio>

#include "utils/chunk_vector.hpp"
#include "utils/allocinfo_allocator.hpp"

template <size_t Size>
struct worker_str
{
  void operator()()
  {
    typedef utils::chunk_vector<std::string, Size> vec_chunk_type;
    typedef std::vector<std::string> vec_type;
    
    vec_chunk_type vec_chunk;
    vec_type vec;

    
    std::cout << "string..." << std::endl;
    
    std::cout << "sizeof: " << sizeof(vec) << " " << sizeof(vec_chunk) << std::endl;

      
    char buffer[1024];
    for (int i = 0; i < 1024 + 100; ++ i) {
      sprintf(buffer, "%06d", i);
      vec.push_back(buffer);
      vec_chunk.push_back(buffer);
    }
    
    std::cerr << "size:" << vec.size() << " " << vec_chunk.size() << std::endl;
    
    for (int i = 0; i < 1024 + 100; ++ i) {
      if (vec[i] != vec_chunk[i])
	std::cerr << "DIFFER!" << std::endl;
      if (*(vec.begin() + i) != *(vec_chunk.begin() + i))
	std::cerr << "DIFFER!" << std::endl;
    }
    
    vec_chunk.clear();
    vec.clear();
    
    std::cerr << "size:" << vec.size() << " " << vec_chunk.size() << std::endl;

    for (int i = 0; i < 1024 + 100; ++ i) {
      sprintf(buffer, "%06d", i);
      vec.push_back(buffer);
      vec_chunk.push_back(buffer);
    }
    
    std::cerr << "size:" << vec.size() << " " << vec_chunk.size() << std::endl;
    
    for (int i = 0; i < 1024 + 100; ++ i) {
      if (vec[i] != vec_chunk[i])
	std::cerr << "DIFFER!" << std::endl;
      if (*(vec.begin() + i) != *(vec_chunk.begin() + i))
	std::cerr << "DIFFER!" << std::endl;
    }

    {
      typename vec_type::iterator iter = vec.begin();
      typename vec_chunk_type::iterator iter_chunk = vec_chunk.begin();
      
      for (/**/; iter_chunk != vec_chunk.end(); ++ iter_chunk, ++ iter) {
	if (*iter != *iter_chunk)
	  std::cerr << "DIFFER!" << std::endl;
      }
    }
    
    std::cerr << "size:" << vec.size() << " " << vec_chunk.size() << std::endl;
    
    for (int i = 0; i < 100; ++ i) {
      typename vec_type::iterator iter = vec.begin();
      typename vec_chunk_type::iterator iter_chunk = vec_chunk.begin();
      
      for (/**/; iter != vec.end(); ++ iter_chunk, ++ iter) {
	if (*iter != *iter_chunk)
	  std::cerr << "DIFFER!" << std::endl;
      }
      
      vec.pop_back();
      vec_chunk.pop_back();
    }
    
    std::cerr << "size:" << vec.size() << " " << vec_chunk.size() << std::endl;
    
    
    vec_chunk.resize(1024 * 4 + 100);
    vec.resize(1024 * 4 + 100);
    
    for (int i = 0; i < 1024 * 4 + 100; ++ i) {
      sprintf(buffer, "%06d", i);
      vec[i] = buffer;
      vec_chunk[i] = buffer;
    }
    
    std::cerr << "size:" << vec.size() << " " << vec_chunk.size() << std::endl;
    
    for (int i = 0; i < 1024 * 4 + 100; ++ i) {
      if (vec[i] != vec_chunk[i])
	std::cerr << "DIFFER!" << std::endl;
      if (*(vec.begin() + i) != *(vec_chunk.begin() + i))
	std::cerr << "DIFFER!" << std::endl;
    }
    
    std::cout << "finished string..." << std::endl;
  }
  
};

template <typename Tp, size_t Size>
struct worker
{
  
  void operator()()
  {
    typedef utils::chunk_vector<Tp, Size, utils::allocinfo_allocator<Tp, std::allocator<Tp> > > vector_type;
    
    utils::allocinfo allocinfo;
  
    std::cout << "worker: initial: " << allocinfo.allocated() << std::endl;

    vector_type vec;
    
    std::cout << "worker: after initialization: " << allocinfo.allocated() << std::endl;
    
    for (int i = 0; i < 1024 + 100; ++ i)
      vec.push_back(i);
    
    std::cout << "worker: after push_back: " << allocinfo.allocated() << std::endl;
    
    std::cout << "sizeof:  " << sizeof(vector_type) << std::endl;
    std::cout << "size: " << vec.size() << std::endl;
    
    {
      typename vector_type::iterator iter = vec.begin();
      for (int i = 0; i < 1024 + 100; ++ i, ++ iter)
	if (*iter != i)
	  std::cerr << "differ?" << std::endl;
    }
    
    {
      int i = 0;
      for (typename vector_type::iterator iter = vec.begin(); iter != vec.end(); ++ iter, ++ i)
	if (*iter != i)
	  std::cerr << "differ?" << std::endl;
      
      if (i != 1024 + 100)
	std::cerr << "differ?" << std::endl;
    }
    
    {
      int vec_size = 1024 + 100;
      for (int i = 0; i < 1024; ++ i) {
	int pos = random() % vec_size;
	if (pos != *(vec.begin() + pos))
	  std::cerr << "differ?" << std::endl;
	
	typename vector_type::iterator iter = vec.begin();
	iter += pos;
	if (pos != *iter)
	  std::cerr << "differ?" << std::endl;
      }
    }
    
    std::cout << "worker: before clear: " << allocinfo.allocated() << std::endl;
    
    vec.clear();
    std::cout << "clear size: " << vec.size() << " " << vec.empty() << std::endl;
    
    std::cout << "worker: terminated: " << allocinfo.allocated() << std::endl;
  }
};

int main(int argc, char** argv)
{
  srandom(time(0));

  utils::allocinfo allocinfo;
  
  std::cout << "initial: " << allocinfo.allocated() << std::endl;

  worker<int, 24>()();
  
  std::cout << "finished int: " << allocinfo.allocated() << std::endl;
  
  worker<uint64_t, 24>()();
  
  std::cout << "finished uint64_t: " << allocinfo.allocated() << std::endl;
  
  worker_str<34>()();
}
