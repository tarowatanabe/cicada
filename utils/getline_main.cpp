//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "utils/getline.hpp"

int main(int argc, char** argv)
{
  std::string line;

  while (utils::getline(std::cin, line))
    std::cout << line << std::endl;

}
