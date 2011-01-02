//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <unistd.h>
#include <time.h>

#include <utils/b_heap.hpp>

int main(int argc, char**argv)
{
  srandom(time(0) * getpid());
  
  utils::b_heap<int> heap;
  std::vector<int> vec;
  
  for (int i = 0; i < 4096; ++ i) {
    vec.push_back(random());
    heap.push(vec.back());
  }
  
  std::sort(vec.begin(), vec.end(), std::greater<int>());

  std::vector<int> vec_heap;
  
  while (! heap.empty()) {
    std::cerr << "value: " << heap.top() << std::endl;
    vec_heap.push_back(heap.top());
    heap.pop();
  }
  

  if (vec != vec_heap)
    std::cerr << "differ" << std::endl;
}
