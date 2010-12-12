
#include <iostream>

#include "alignment.hpp"

int main(int argc, char** argv)
{
  typedef cicada::Alignment alignment_type;
  
  alignment_type alignment("5-6 7-8 0-2 5-6");
  for (alignment_type::const_iterator iter = alignment.begin(); iter != alignment.end(); ++ iter)
    std::cout << "align: " << *iter << std::endl;
  std::cout << alignment << std::endl;

  alignment_type::point_type point("500-40");
  std::cout << "point: " << point << std::endl;

  alignment_type::point_type point2("500- 40");
  std::cout << "point: " << point2 << std::endl;
}
