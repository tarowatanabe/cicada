//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/prefix.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/bytestream.h>

namespace cicada
{
  namespace stemmer
  {
    Stemmer::symbol_type Prefix::operator[](const symbol_type& word) const
    {
      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not prefixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;

      symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
    
      if (word.id() >= __cache.size())
	__cache.resize(word.id() + 1, vocab_type::EMPTY);

      if (__cache[word.id()] == vocab_type::EMPTY) {
	UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
    
	const size_t index = uword.moveIndex32(0, static_cast<int>(size));
      
	UnicodeString uword_prefix;
	uword.extractBetween(0, index, uword_prefix);
      
	if (uword_prefix.length() < uword.length()) {
	  uword_prefix.append('+');
	
	  std::string word_prefix;
	  StringByteSink<std::string> __sink(&word_prefix);
	  uword_prefix.toUTF8(__sink);
	
	  __cache[word.id()] = word_prefix;
	} else
	  __cache[word.id()] = word;
      }
    
      return __cache[word.id()];
    }

  };
};
