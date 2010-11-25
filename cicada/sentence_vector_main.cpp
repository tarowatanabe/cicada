//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "sentence_vector.hpp"

int main(int argc, char** argv)
{
  typedef cicada::Sentence sentence_type;
  typedef cicada::SentenceVector sentence_set_type;

  sentence_set_type sentences("Good morning |||  bad evening |||");
  
  sentence_set_type::const_iterator siter_end = sentences.end();
  for (sentence_set_type::const_iterator siter = sentences.begin(); siter != siter_end; ++ siter)
    std::cout << "sentence: " <<  *siter << ':' << '\n';
  std::cout << "END" << '\n';
  std::cout << "sentences: " << sentences << ':' << std::endl;
}



