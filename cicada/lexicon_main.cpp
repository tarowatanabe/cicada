//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "lexicon.hpp"
#include "sentence.hpp"

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << argv[0] << " lexicon-file" << std::endl;
    return 1;
  }
  
  cicada::Lexicon lexicon(argv[1]);
  
  cicada::Sentence sentence;
  while (std::cin >> sentence)
    if (! sentence.empty())
      std::cerr << "value: " << lexicon(sentence.begin(), sentence.end() - 1, sentence.back()) << std::endl
		<< "exists: " << lexicon.exists(sentence.begin(), sentence.end() - 1) << std::endl
		<< "exists2: " << lexicon.exists(sentence.begin(), sentence.end() - 1, sentence.back()) << std::endl;
  
}
