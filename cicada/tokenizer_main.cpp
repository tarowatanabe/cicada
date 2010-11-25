//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "tokenizer.hpp"

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << argv[0] << " tokenizer-spec" << std::endl;
    std::cout << cicada::Tokenizer::lists();
    return 1;
  }
  
  cicada::Tokenizer& tokenizer(cicada::Tokenizer::create(argv[1]));
  
  cicada::Sentence sentence;
  cicada::Sentence tokenized;
  while (std::cin >> sentence) {
    std::cout << "original: " << sentence << std::endl;
    
    tokenizer(sentence, tokenized);
    
    std::cout << "tokenized: " << tokenized << std::endl;
  }
}
