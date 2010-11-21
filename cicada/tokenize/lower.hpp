// -*- mode: c++ -*-

#ifndef __CICADA__TOKENIZE__LOWER__HPP__
#define __CICADA__TOKENIZE__LOWER__HPP__ 1

#include <vector>
#include <string>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>
#include <unicode/bytestream.h>

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace tokenize
  {
    
#ifndef HAVE_TLS
    static void __lower_no_clean_up_func(cicada::Stemmer* x)
    {
      
    }
#endif
    
    template <typename Sent, typename Tokenized>
    void lower(const Sent& sentence, Tokenized& __tokenized)
    {
      typedef Sent sentence_type;
      typedef typename sentence_type::value_type word_type;
      
      typedef cicada::Stemmer stemmer_type;

#ifdef HAVE_TLS
      static __thread stemmer_type* __stemmer_tls = 0;
      if (! __stemmer_tls)
	__stemmer_tls = &stemmer_type::create("lower");
      stemmer_type& stemmer = *__stemmer_tls;
#else
      static boost::thread_specific_ptr<stemmer_type> __stemmer(__lower_no_clean_up_func);
      if (! __stemmer.get())
	__stemmer.reset(&stemmer_type::create("lower"));
      
      stemmer_type& stemmer = *__stemmer;
#endif
      
      std::vector<word_type, std::allocator<word_type> > tokenized;
      
      typename sentence_type::const_iterator siter_end = sentence.end();
      for (typename sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter)
	tokenized.push_back(stemmer(*siter));
      
      __tokenized = Tokenized(tokenized.begin(), tokenized.end());
    }
    
  };
};

#endif
