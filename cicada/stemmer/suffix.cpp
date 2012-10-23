//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/suffix.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>

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


      symbol_pair_set_type& __cache = const_cast<symbol_pair_set_type&>(cache);
      symbol_pair_type& pair = __cache[word.id() & (__cache.size() - 1)];
    
      if (pair.first != word) {
	icu::UnicodeString uword = icu::UnicodeString::fromUTF8(static_cast<const std::string&>(word));
      
	const size_t index = uword.moveIndex32(uword.length(), - int(size));
      
	icu::UnicodeString uword_suffix;
	uword.extractBetween(index, uword.length(), uword_suffix);
      
	if (uword_suffix.length () < uword.length()) {
	  uword_suffix.insert(0, '+');
	
	  std::string word_suffix;
	  uword_suffix.toUTF8String(word_suffix);
	
	  pair.second = word_suffix;
	} else
	  pair.second = word;
	
	pair.first = word;
      }
    
      return pair.second;
    }

  };
};
