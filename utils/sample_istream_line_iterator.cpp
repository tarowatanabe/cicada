//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <utils/istream_line_iterator.hpp>

int main(int argc, char** argv)
{
  
  size_t num_lines = 0;
  utils::istream_line_iterator iter_begin(std::cin);
  utils::istream_line_iterator iter_end;
  for (utils::istream_line_iterator iter = iter_begin; iter != iter_end; ++ iter)
    ++ num_lines;
  
  std::cout << "# of lines: " << num_lines << std::endl;
}
