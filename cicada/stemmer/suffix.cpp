//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/suffix.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>

namespace cicada
{
  namespace stemmer
  {
    std::string Suffix::operator()(const utils::piece& word) const
    {
      if (word.empty()) return word;
      
      const size_type word_size = word.size();
      
      // SGML-like symbols are not suffixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      icu::UnicodeString uword = icu::UnicodeString::fromUTF8(icu::StringPiece(word.data(), word.size()));
      
      const size_t index = uword.moveIndex32(uword.length(), - int(size));
      
      icu::UnicodeString uword_suffix;
      uword.extractBetween(index, uword.length(), uword_suffix);
      
      if (uword_suffix.length () < uword.length()) {
	uword_suffix.insert(0, '+');
	
	std::string word_suffix;
	uword_suffix.toUTF8String(word_suffix);
	
	return word_suffix;
      } else
	return word;
    }

  };
};
