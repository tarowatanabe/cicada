//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <iostream>

#include "matcher.hpp"

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << argv[0] << " matcher-spec" << std::endl;
    std::cout << cicada::Matcher::lists();
    return 1;
  }

  try {
    cicada::Matcher& matcher(cicada::Matcher::create(argv[1]));
    
    std::string word1;
    std::string word2;
    while (std::cin >> word1 >> word2)
      std::cout << "word1: " << word1 << " word2: " << word2 << " matched? " << matcher(word1, word2) << std::endl;
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}
