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
      std::auto_ptr<icu::Transliterator> trans(icu::Transliterator::createInstance(icu::UnicodeString::fromUTF8("Lower"), UTRANS_FORWARD, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      
      pimpl = trans.release();
    }
    
    English::~English()
    {
      std::auto_ptr<icu::Transliterator> tmp(static_cast<icu::Transliterator*>(pimpl));
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
	icu::UnicodeString uword = icu::UnicodeString::fromUTF8(static_cast<const std::string&>(word));

	std::string signature = "<UNK";
	
	// signature for English, taken from Stanford parser's getSignature5
	int num_caps = 0;
	bool has_digit = false;
	bool has_dash = false;
	bool has_lower = false;
	
	size_t length = 0;
	UChar32 ch0 = 0;
	UChar32 ch_1 = 0;
	UChar32 ch_2 = 0;
	icu::StringCharacterIterator iter(uword);
	for (iter.setToStart(); iter.hasNext(); ++ length) {
	  UChar32 ch = iter.next32PostInc();
	  
	  // keep initial char...
	  if (ch0 == 0)
	    ch0 = ch;
	  
	  ch_2 = ch_1;
	  ch_1 = ch;
	  
	  if (u_getNumericValue(ch) != U_NO_NUMERIC_VALUE)
	    has_digit = true;
	  if (u_hasBinaryProperty(ch,  UCHAR_DASH))
	    has_dash = true;
	  if (u_isUAlphabetic(ch)) {
	    if (u_isULowercase(ch))
	      has_lower = true;
	    else if (u_istitle(ch)) {
	      has_lower = true;
	      ++num_caps;
	    } else
	      ++num_caps;
	  }
	}
	
	// transform into lower...
	static_cast<icu::Transliterator*>(pimpl)->transliterate(uword);
	ch_2 = (ch_2 ? u_tolower(ch_2) : ch_2);
	ch_1 = (ch_1 ? u_tolower(ch_1) : ch_1);
	
	// we do not check loc...
	if (u_isUUppercase(ch0) || u_istitle(ch0))
	  signature += "-CAPS";
	else if (! u_isUAlphabetic(ch0) && num_caps > 0)
	  signature += "-CAPS";
	else if (has_lower)
	  signature += "-LC";
	
	if (has_digit)
	  signature += "-NUM";
	if (has_dash)
	  signature += "-DASH";
	
	if (length >= 3 && ch_1 == 's') {
	  if (ch_2 != 's' && ch_2 != 'i' && ch_2 != 'u')
	    signature += "-s";
	} else if (length >= 5 && ! has_dash && ! (has_digit && num_caps > 0)) {
	  if (uword.endsWith("ed"))
	    signature += "-ed";
	  else if (uword.endsWith("ing"))
	    signature += "-ing";
	  else if (uword.endsWith("ion"))
	    signature += "-ion";
	  else if (uword.endsWith("er"))
	    signature += "-er";
	  else if (uword.endsWith("est"))
	    signature += "-est";
	  else if (uword.endsWith("ly"))
	    signature += "-ly";
	  else if (uword.endsWith("ity"))
	    signature += "-ity";
	  else if (uword.endsWith("y"))
	    signature += "-y";
	  else if (uword.endsWith("al"))
	    signature += "-al";
	}

	signature += '>';
	
	__cache[word.id()] = signature;
      }
      
      return __cache[word.id()];
    }
      
  };
};
