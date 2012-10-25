//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/prefix.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>

namespace cicada
{
  namespace stemmer
  {
    std::string Prefix::operator()(const utils::piece& word) const
    {
      if (word.empty()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not prefixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      icu::UnicodeString uword = icu::UnicodeString::fromUTF8(icu::StringPiece(word.data(), word.size()));
      
      const size_t index = uword.moveIndex32(0, static_cast<int>(size));
      
      icu::UnicodeString uword_prefix;
      uword.extractBetween(0, index, uword_prefix);
      
      if (uword_prefix.length() < uword.length()) {
	uword_prefix.append('+');
	
	std::string word_prefix;
	uword_prefix.toUTF8String(word_prefix);
	
	return word_prefix;
      } else
	return word;
    }

  };
};
