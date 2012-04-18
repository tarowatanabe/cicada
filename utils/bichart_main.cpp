//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//
#include <cstdlib>
#include <algorithm>
#include <iostream>

#include "utils/bichart.hpp"

int main(int argc, char** argv)
{
  utils::bichart<int> chart;
  
  srandom(time(0) * getpid());
  
  for (int j = 0; j < 10; ++ j) {
    int dim1 = 0;
    int dim2 = 0;
    while (dim1 == 0)
      dim1 = random() % 100;

    while (dim2 == 0)
      dim2 = random() % 100;
    
    std::cout << "dim1=" << dim1 << " dim2="<< dim2 << std::endl;
    
    chart.resize(dim1, dim2);

    std::vector<int*> pointers;
    
    for (int first1 = 0; first1 != dim1; ++ first1)
      for (int last1 = first1; last1 != dim1; ++ last1)
	for (int first2 = 0; first2 != dim2; ++ first2)
	  for (int last2 = first2; last2 != dim2; ++ last2) {
	    chart(first1, last1, first2, last2) = random();
	    
	    pointers.push_back(&chart(first1, last1, first2, last2));
	  }
    
    std::sort(pointers.begin(), pointers.end());
    for (size_t pos = 1; pos != pointers.size(); ++ pos)
      if (pointers[pos - 1] == pointers[pos])
	std::cerr << "duplicated position?" << std::endl;
  }
}
