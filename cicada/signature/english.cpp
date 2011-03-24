//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "signature/english.hpp"

#include <boost/thread.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>
#include <unicode/schriter.h>

#include "utils/config.hpp"

namespace cicada
{
  namespace signature
  {
    English::English() : pimpl(0) 
    {
      UErrorCode status = U_ZERO_ERROR;
      std::auto_ptr<Transliterator> trans(Transliterator::createInstance(UnicodeString::fromUTF8("Lower"), UTRANS_FORWARD, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      
      pimpl = trans.release();
    }
    
    English::~English()
    {
      std::auto_ptr<Transliterator> tmp(static_cast<Transliterator*>(pimpl));
    }
    
    Signature::symbol_type English::operator[](const symbol_type& word) const
    {
      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
      
      const size_type word_size = word.size();
    
      // SGML-like symbols are used as-is
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
      
      if (word.id() >= __cache.size())
	__cache.resize(word.id() + 1, vocab_type::EMPTY);
    
      if (__cache[word.id()] == vocab_type::EMPTY) {
	UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
	
	// signature for English, taken from Stanford parser's getSignature5
	int num_caps = 0;
	bool has_digit = false;
	bool has_dash = false;
	bool has_lower = false;
	
	StringCharacterIterator iter(uword);
	for (iter.setToStart(); iter.hasNext(); /**/) {
	  UChar32 ch = iter.next32PostInc();
	  
	}
	
	// transform into lower...
	static_cast<Transliterator*>(pimpl)->transliterate(uword);
      }
      
      return __cache[word.id()];
    }
      
  };
};
