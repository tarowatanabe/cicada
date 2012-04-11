//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//
#include <cstdlib>

#include <iostream>

#include "utils/chart.hpp"

int main(int argc, char** argv)
{
  utils::chart<int> chart;
  
  srandom(time(0) * getpid());
  
  for (int j = 0; j < 10; ++ j) {
    int dim1 = 0;
    while (dim1 == 0)
      dim1 = random() % 100;
    
    std::cout << dim1 << std::endl;
  
    chart.resize(dim1);

    std::vector<int*> pointers;
    
    for (int first = 0; first != dim1; ++ first)
      for (int last = first; last != dim1; ++ last)
	pointers.push_back(&chart(first, last));
    
    std::sort(pointers.begin(), pointers.end());
    for (size_t pos = 1; pos != pointers.size(); ++ pos)
      if (pointers[pos - 1] == pointers[pos])
	std::cerr << "duplicated position?" << std::endl;
  }
}
