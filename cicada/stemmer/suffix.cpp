//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/suffix.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/bytestream.h>

namespace cicada
{
  namespace stemmer
  {
    Stemmer::symbol_type Suffix::operator[](const symbol_type& word) const
    {
      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not suffixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;

      symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
    
      if (word.id() >= __cache.size())
	__cache.resize(word.id() + 1, vocab_type::EMPTY);
    
      if (__cache[word.id()] == vocab_type::EMPTY) {
	UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
      
	const size_t index = uword.moveIndex32(uword.length(), - int(size));
      
	UnicodeString uword_suffix;
	uword.extractBetween(index, uword.length(), uword_suffix);
      
	if (uword_suffix.length () < uword.length()) {
	  uword_suffix.insert(0, '+');
	
	  std::string word_suffix;
	  StringByteSink<std::string> __sink(&word_suffix);
	  uword_suffix.toUTF8(__sink);
	
	  __cache[word.id()] = word_suffix;
	} else
	  __cache[word.id()] = word;
      }
    
      return __cache[word.id()];
    }

  };
};
