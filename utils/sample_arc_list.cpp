//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdlib.h>

#include <iostream>
#include <vector>

#include "utils/arc_list.hpp"
#include "utils/atomicop.hpp"

#include <boost/thread.hpp>

typedef utils::arc_list<int, int, 32> cache_type;

int main(int argc, char** argv)
{
  cache_type cache;
  
  for (int i = 0; i < 1024; ++ i) {
    std::pair<cache_type::iterator, bool> result = cache.find(i);
    
    if (! result.second)
      result.first->second = i;

    if (result.first->first != i)
      std::cerr << "DIFFERENT KEY? " << result.first->first << " " << i << std::endl;
    
    std::cerr << "result: " << result.first->first
	      << " " << result.first->second
	      << " " << result.second
	      << std::endl;
  }
  
  std::cerr << "RANDOM!" << std::endl;
  
  srandom(getpid() * time(0));
  
  for (int j = 0; j < 1024 * 5; ++ j) {
    int i = random() % 64;
    
    std::pair<cache_type::iterator, bool> result = cache.find(i);
    
    if (! result.second)
      result.first->second = i;
    else if (result.first->second != i)
      std::cerr << "DIFFERENT cache value???" << std::endl;
    
    if (result.first->first != i)
      std::cerr << "DIFFERENT KEY? " << result.first->first << " " << i << std::endl;
    
    std::cerr << "result: " << result.first->first
	      << " " << result.first->second
	      << " " << result.second
	      << std::endl;
  }
  

}
