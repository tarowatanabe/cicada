// -*- mode: c++ -*-

#ifndef __CICADA__TOKENIZER__STEMMER__HPP__
#define __CICADA__TOKENIZER__STEMMER__HPP__ 1

#include <vector>
#include <string>

#include <cicada/stemmer.hpp>
#include <cicada/tokenizer.hpp>

namespace cicada
{
  namespace tokenizer
  {
    class Stemmer : public cicada::Tokenizer
    {
    public:
      Stemmer(const cicada::Stemmer* __stemmer) : stemmer(__stemmer) {}
      
    protected:
      virtual void tokenize(const sentence_type& source, sentence_type& tokenized) const
      {
	tokenized.clear();
	sentence_type::const_iterator siter_end = source.end();
	for (sentence_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter)
	  tokenized.push_back(stemmer->operator()(*siter));
      }
      
    private:
      const cicada::Stemmer* stemmer;
    };
  };
};

#endif
