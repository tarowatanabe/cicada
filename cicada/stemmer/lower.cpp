//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/lower.hpp"

#include <boost/thread.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>
#include <unicode/bytestream.h>

#include "utils/config.hpp"

namespace cicada
{
  namespace stemmer
  {
    Lower::Lower() : pimpl(0)
    {
      UErrorCode status = U_ZERO_ERROR;
      std::auto_ptr<Transliterator> trans(Transliterator::createInstance(UnicodeString::fromUTF8("Lower"), UTRANS_FORWARD, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      
      pimpl = trans.release();
    }
    
    Lower::~Lower()
    {
      std::auto_ptr<Transliterator> tmp(static_cast<Transliterator*>(pimpl));
    }
    
    Stemmer::symbol_type Lower::operator[](const symbol_type& word) const
    {
      if (! pimpl)
	throw std::runtime_error("no lower caser?");

      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not lowered...
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
      
      if (word.id() >= __cache.size())
	__cache.resize(word.id() + 1, vocab_type::EMPTY);
    
      if (__cache[word.id()] == vocab_type::EMPTY) {
	std::string word_lower;
	UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
	
	static_cast<Transliterator*>(pimpl)->transliterate(uword);
	
	StringByteSink<std::string> __sink(&word_lower);
	uword.toUTF8(__sink);
	
	__cache[word.id()] = word_lower;
      }
    
      return __cache[word.id()];
    }

  };
};
