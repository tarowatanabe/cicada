//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/prefix.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>

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
      
      symbol_pair_set_type& __cache = const_cast<symbol_pair_set_type&>(cache);
      symbol_pair_type& pair = __cache[word.id() & (__cache.size() - 1)];
      
      if (pair.first != word) {
	icu::UnicodeString uword = icu::UnicodeString::fromUTF8(static_cast<const std::string&>(word));
	
	const size_t index = uword.moveIndex32(0, static_cast<int>(size));
	
	icu::UnicodeString uword_prefix;
	uword.extractBetween(0, index, uword_prefix);
	
	if (uword_prefix.length() < uword.length()) {
	  uword_prefix.append('+');
	  
	  std::string word_prefix;
	  uword_prefix.toUTF8String(word_prefix);
	  
	  pair.second = word_prefix;
	} else
	  pair.second = word;
	
	pair.first = word;
      }
      
      return pair.second;
    }

  };
};
