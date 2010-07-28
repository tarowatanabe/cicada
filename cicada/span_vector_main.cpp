
#include <iostream>

#include "span_vector.hpp"

int main(int argc, char** argv)
{
  cicada::SpanVector input;
  while (std::cin >> input)
    std::cout << input << std::endl;
}
