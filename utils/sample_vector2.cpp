//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//
#include <cstdlib>

#include <iostream>

#include "utils/vector2.hpp"
#include "utils/vector3.hpp"

#include "utils/vector2_aligned.hpp"
#include "utils/vector3_aligned.hpp"

int main(int argc, char** argv)
{
  utils::vector2<int> vec2;
  utils::vector3<int> vec3;

  utils::vector2_aligned<int> vec2_aligned;
  utils::vector3_aligned<int> vec3_aligned;

  srandom(time(0) * getpid());
  
  for (int j = 0; j < 10; ++ j) {
  
    int dim1 = 0;
    int dim2 = 0;
    int dim3 = 0;
  
    while (dim1 == 0 || dim2 == 0 || dim3 == 0) {
      dim1 = random() % 100;
      dim2 = random() % 100;
      dim3 = random() % 100;
    }
    
    std::cout << dim1 << ' ' << dim2 << ' ' << dim3 << std::endl;
  
    vec2.resize(dim1, dim2);
    vec3.resize(dim1, dim2, dim3);

    vec2_aligned.resize(dim1, dim2);
    vec3_aligned.resize(dim1, dim2, dim3);
  
    for (int i = 0; i < 1000; ++ i) {
      const int pos1 = random() % dim1;
      const int pos2 = random() % dim2;
      const int pos3 = random() % dim3;
    
      if (&(vec2(pos1, 0)) != &(*vec2.begin(pos1)))
	std::cerr << "DIFFER..." << std::endl;
      if (&(vec2(pos1, pos2)) != &(*(vec2.begin(pos1) + pos2)))
	std::cerr << "DIFFER..." << std::endl;
      if (&(vec2(pos1, pos2)) != &(vec2[pos1][pos2]))
	std::cerr << "DIFFER..." << std::endl;
      
      if (&(vec2_aligned(pos1, 0)) != &(*vec2_aligned.begin(pos1)))
	std::cerr << "DIFFER..." << std::endl;
      if (&(vec2_aligned(pos1, pos2)) != &(*(vec2_aligned.begin(pos1) + pos2)))
	std::cerr << "DIFFER..." << std::endl;
      if (&(vec2_aligned(pos1, pos2)) != &(vec2_aligned[pos1][pos2]))
	std::cerr << "DIFFER..." << std::endl;
      
      for (int pos = 0; pos < dim2; ++ pos) {
	if (&(vec2(pos1, pos)) != &(vec2[pos1][pos]))
	  std::cerr << "DIFFER..." << std::endl;
	
	if (&(vec2_aligned(pos1, pos)) != &(vec2_aligned[pos1][pos]))
	  std::cerr << "DIFFER..." << std::endl;
      }
      
      if (&(vec3(pos1, 0, 0)) != &(*vec3.begin(pos1)))
	std::cerr << "DIFFER..." << std::endl;
      if (&(vec3(pos1, pos2, 0)) != &(*(vec3.begin(pos1, pos2))))
	std::cerr << "DIFFER..." << std::endl;
      if (&(vec3(pos1, pos2, pos3)) != &(*(vec3.begin(pos1, pos2) + pos3)))
	std::cerr << "DIFFER..." << std::endl;
      
      if (&(vec3_aligned(pos1, 0, 0)) != &(*vec3_aligned.begin(pos1)))
	std::cerr << "DIFFER..." << std::endl;
      if (&(vec3_aligned(pos1, pos2, 0)) != &(*(vec3_aligned.begin(pos1, pos2))))
	std::cerr << "DIFFER..." << std::endl;
      if (&(vec3_aligned(pos1, pos2, pos3)) != &(*(vec3_aligned.begin(pos1, pos2) + pos3)))
	std::cerr << "DIFFER..." << std::endl;
    }
  }
}
