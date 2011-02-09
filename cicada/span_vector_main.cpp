//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "span_vector.hpp"

int main(int argc, char** argv)
{
  std::cout << "span: " << cicada::SpanVector::Span("1-2") << std::endl;
  std::cout << "span: " << cicada::SpanVector::Span("1-1:good") << std::endl;

  cicada::SpanVector input;
  while (std::cin >> input) {
    std::cout << "input: " << input << std::endl;
    
    for (size_t i = 0; i != input.size(); ++ i)
      std::cout << "span: " << input[i] << std::endl;
  }
}
