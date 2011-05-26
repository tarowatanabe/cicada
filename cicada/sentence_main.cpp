//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "sentence.hpp"

int main(int argc, char** argv)
{
  typedef cicada::Sentence sentence_type;
  
  sentence_type sentence("good morning boy");
  
  sentence_type::const_iterator siter_end = sentence.end();
  for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter)
    std::cout << "word: " <<  *siter << '\n';
  std::cout << "END" << '\n';

  std::cout << "generation: " << sentence << std::endl;
  
  std::string parsing("good morning boy , yes ||| bood");
  std::string::const_iterator iter = parsing.begin();
  std::string::const_iterator end = parsing.end();

  if (sentence.assign(iter, end)) {
    sentence_type::const_iterator siter_end = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter)
      std::cout << "word: " <<  *siter << '\n';
    std::cout << "END" << '\n';

    std::cout << "remain: " << std::string(iter, end) << std::endl;
    
  } else
    std::cout << "failed:" << std::endl;
  
  
  std::cout << "generation: " << sentence << std::endl;
  

  while (std::cin >> sentence)
    std::cout << sentence << std::endl;
}
