//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//


#include "format.hpp"

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << argv[0] << " format-spec" << std::endl;
    std::cout << cicada::Format::lists();
    return 1;
  }
  
  cicada::Format& format(cicada::Format::create(argv[1]));
  
  std::string line;
  cicada::Format::phrase_set_type results;
  while (std::getline(std::cin, line)) {
    format(line, results);
    
    for (size_t i = 0; i != results.size(); ++ i)
      std::cout << i << ": " << results[i] << std::endl;
  }
}
