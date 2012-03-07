//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <queue>

#include <unistd.h>
#include <time.h>

#include <utils/resource.hpp>
#include <utils/b_heap.hpp>
#include <utils/std_heap.hpp>

template <typename Heap>
void process()
{
  Heap heap;

  heap.reserve(1024 * 1024);

  for (int i = 0; i < 1024 * 1024; ++ i)
    heap.push(random());

  for (int i = 0; i < 1024 * 1024; ++ i)
    if (heap.empty())
      heap.push(random());
    else if (random() & 1) {
      heap.top();
      heap.pop();
    } else
      heap.push(random());

  while (! heap.empty()) {
    const int top = heap.top();
    heap.pop();
  }
}

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
    //std::cerr << "value: " << heap.top() << std::endl;
    vec_heap.push_back(heap.top());
    heap.pop();
  }
  

  if (vec != vec_heap)
    std::cerr << "differ" << std::endl;

  utils::resource std_heap_start;
  process<utils::std_heap<int> >();
  utils::resource std_heap_end;

  utils::resource b_heap_start;
  process<utils::b_heap<int, std::vector<int>, std::less<int>, 512 /sizeof(int)> >();
  utils::resource b_heap_end;
  
  
    std::cout << "b heap cpu time: " << (b_heap_end.cpu_time() - b_heap_start.cpu_time())
	      << " user time: " << (b_heap_end.user_time() - b_heap_start.user_time())
	      << std::endl;
    std::cout << "std heap cpu time: " << (std_heap_end.cpu_time() - std_heap_start.cpu_time())
	      << " user time: " << (std_heap_end.user_time() - std_heap_start.user_time())
	      << std::endl;
}
