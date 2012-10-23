//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/simplified.hpp"

#include <boost/thread.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>

#include "utils/config.hpp"

namespace cicada
{
  namespace stemmer
  {
    Simplified::Simplified() : pimpl(0)
    {
      UErrorCode status = U_ZERO_ERROR;
      std::auto_ptr<icu::Transliterator> trans(icu::Transliterator::createInstance(icu::UnicodeString::fromUTF8("Traditional-Simplified"), UTRANS_FORWARD, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      
      pimpl = trans.release();
    }
    
    Simplified::~Simplified()
    {
      std::auto_ptr<icu::Transliterator> tmp(static_cast<icu::Transliterator*>(pimpl));
    }
    
    Stemmer::symbol_type Simplified::operator[](const symbol_type& word) const
    {
      if (! pimpl)
	throw std::runtime_error("no simplified?");

      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not simplifieded...
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      symbol_pair_set_type& __cache = const_cast<symbol_pair_set_type&>(cache);
      symbol_pair_type& pair = __cache[word.id() & (__cache.size() - 1)];
    
      if (pair.first != word) {
	icu::UnicodeString uword = icu::UnicodeString::fromUTF8(static_cast<const std::string&>(word));
	
	static_cast<icu::Transliterator*>(pimpl)->transliterate(uword);
	
	std::string word_simplified;
	uword.toUTF8String(word_simplified);
	
	pair.first  = word;
	pair.second = word_simplified;
      }
    
      return pair.second;
    }

  };
};
