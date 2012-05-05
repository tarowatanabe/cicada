//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/fullwidth.hpp"

#include <boost/thread.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>

#include "utils/config.hpp"

namespace cicada
{
  namespace stemmer
  {
    Fullwidth::Fullwidth() : pimpl(0)
    {
      UErrorCode status = U_ZERO_ERROR;
      std::auto_ptr<icu::Transliterator> trans(icu::Transliterator::createInstance(icu::UnicodeString::fromUTF8("Halfwidth-Fullwidth"), UTRANS_FORWARD, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      
      pimpl = trans.release();
    }
    
    Fullwidth::~Fullwidth()
    {
      std::auto_ptr<icu::Transliterator> tmp(static_cast<icu::Transliterator*>(pimpl));
    }
    
    Stemmer::symbol_type Fullwidth::operator[](const symbol_type& word) const
    {
      if (! pimpl)
	throw std::runtime_error("no fullwidth?");

      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not fullwidthed...
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
      
      if (word.id() >= __cache.size())
	__cache.resize(word.id() + 1, vocab_type::EMPTY);
    
      if (__cache[word.id()] == vocab_type::EMPTY) {
	
	icu::UnicodeString uword = icu::UnicodeString::fromUTF8(static_cast<const std::string&>(word));
	
	static_cast<icu::Transliterator*>(pimpl)->transliterate(uword);
	
	std::string word_fullwidth;
	uword.toUTF8String(word_fullwidth);
	
	__cache[word.id()] = word_fullwidth;
      }
    
      return __cache[word.id()];
    }

  };
};
