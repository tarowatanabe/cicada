// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TOKENIZER__TOKENIZE__HPP__
#define __CICADA__TOKENIZER__TOKENIZE__HPP__ 1

#include <vector>

#include <cicada/tokenizer.hpp>

namespace cicada
{
  namespace tokenizer
  {
    class Tokenize : public cicada::Tokenizer
    {
    private:
      typedef cicada::Tokenizer tokenizer_type;
      typedef std::vector<const tokenizer_type*, std::allocator<const tokenizer_type*> > tokenizer_set_type;

    public:
      void insert(const Tokenizer& tokenizer)
      {
	tokenizers.push_back(&tokenizer);
      }

      bool empty() const { return tokenizers.empty(); }
      
    protected:
      virtual void tokenize(const sentence_type& source, sentence_type& tokenized) const
      {
	sentence_type sentence;
	tokenized = source;
	tokenizer_set_type::const_iterator titer_end = tokenizers.end();
	for (tokenizer_set_type::const_iterator titer = tokenizers.begin(); titer != titer_end; ++ titer) {
	  sentence.swap(tokenized);
	  (*titer)->operator()(sentence, tokenized);
	}
      }
      
    private:
      tokenizer_set_type tokenizers;
    };
  };
};

#endif
