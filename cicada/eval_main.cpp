//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "eval.hpp"
#include "sentence_vector.hpp"

typedef cicada::eval::Scorer   scorer_type;
typedef cicada::Sentence sentence_type;
typedef cicada::SentenceVector sentence_set_type;

void process(const std::string& name, const sentence_set_type& refset, const sentence_type& tstset)
{
  scorer_type::scorer_ptr_type scorer(scorer_type::create(name));
  scorer->insert(refset);
  
  scorer_type::score_ptr_type score = scorer->score(tstset);

  std::cout << "score: " <<score->score() << " desc: " << *score << std::endl;
  
  const std::string encoded = score->encode();

  std::cout << "encoded: " << encoded << std::endl;

  scorer_type::score_ptr_type decoded = scorer_type::score_type::decode(encoded);
  
  if (! decoded)
    std::cerr << "decoding failed" << std::endl;
  else if (*decoded != *score)
    std::cerr << "encoding/decoding failed" << std::endl;
}

int main(int argc, char** argv)
{
  sentence_set_type refset("aaa bbb ccc ddd ||| aaa ccc");
  sentence_type     tstset("aaa ccc ddd");

  process("bleu", refset, tstset);
  process("ter", refset, tstset);
  process("wer", refset, tstset);
  process("per", refset, tstset);
  process("sk", refset, tstset);
  process("sb", refset, tstset);
  process("wlcs", refset, tstset);
  process("combined:reward=true,metric=bleu,weight=0.5,metric=ter,weight=-0.5", refset, tstset);
}
