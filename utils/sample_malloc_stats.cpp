//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>

#include "malloc_stats.hpp"

int main(int argc, char** argv)
{
  typedef std::vector<int> vec_type;
  typedef std::vector<vec_type> vecs_type;

  vecs_type vecs;
  for (int i = 0; i < 1024 * 4; ++ i) {
    vecs.push_back(vec_type(1024));
    
    std::cout << "allocated: " << utils::malloc_stats::allocated()
	      << " used: " << utils::malloc_stats::used()
	      << std::endl;
  }

  vecs.clear();
  
  std::cout << "allocated: " << utils::malloc_stats::allocated()
	    << " used: " << utils::malloc_stats::used()
	    << std::endl;
  
  vecs_type(vecs).swap(vecs);
  
  std::cout << "allocated: " << utils::malloc_stats::allocated()
	    << " used: " << utils::malloc_stats::used()
	    << std::endl;
  
}
