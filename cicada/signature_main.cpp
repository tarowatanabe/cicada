//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "signature.hpp"

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << argv[0] << " signature-spec" << std::endl;
    std::cout << cicada::Signature::lists();
    return 1;
  }

  
  
  cicada::Signature& signature(cicada::Signature::create(argv[1]));
  
  std::string word;
  while (std::cin >> word)
    std::cout << "word: " << word << " signature: " << signature[word] << std::endl;

  
}
