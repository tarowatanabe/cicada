//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "alignment.hpp"

#ifdef HAVE_MSGPACK_HPP
#include <msgpack.hpp>
#include <cicada/msgpack/alignment.hpp>
#endif

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
