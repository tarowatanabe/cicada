// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TOKENIZER__LOWER__HPP__
#define __CICADA__TOKENIZER__LOWER__HPP__ 1

#include <vector>
#include <string>

#include <cicada/stemmer.hpp>
#include <cicada/tokenizer.hpp>

namespace cicada
{
  namespace tokenizer
  {
    class Lower : public cicada::Tokenizer
    {
    public:
      Lower() : lower(&cicada::Stemmer::create("lower")) {}
      
    protected:
      virtual void tokenize(const sentence_type& source, sentence_type& tokenized) const
      {
	tokenized.clear();
	sentence_type::const_iterator siter_end = source.end();
	for (sentence_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter)
	  tokenized.push_back(lower->operator()(*siter));
      }
      
    private:
      cicada::Stemmer* lower;
    };
  };
};

#endif
