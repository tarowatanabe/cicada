//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/digit.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>

namespace cicada
{
  namespace stemmer
  {
    Stemmer::symbol_type Digit::operator[](const symbol_type& word) const
    {
      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not prefixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
    
      symbol_pair_set_type& __cache = const_cast<symbol_pair_set_type&>(cache);
      symbol_pair_type& pair = __cache[word.id() & (__cache.size() - 1)];
    
      if (pair.first != word) {
	pair.first = word;

	icu::UnicodeString uword = icu::UnicodeString::fromUTF8(static_cast<const std::string&>(word));
      
	if (u_getIntPropertyValue(uword.char32At(0), UCHAR_NUMERIC_TYPE) == U_NT_NONE
	    && u_getIntPropertyValue(uword.char32At(uword.length() - 1), UCHAR_NUMERIC_TYPE) == U_NT_NONE)
	  pair.second = word;
	else {
	  bool found = false;
	  icu::UnicodeString uword_digits("<digit-");
	  icu::StringCharacterIterator iter(uword);
	  for (iter.setToStart(); iter.hasNext(); /**/) {
	    UChar32 c = iter.next32PostInc();
	  
	    //const int32_t numeric_type = u_getIntPropertyValue(c, UCHAR_NUMERIC_TYPE);
	    const double numeric_value = u_getNumericValue(c);
	    const int32_t numeric_int = int(numeric_value);
	  
	    const bool replace =(numeric_value != U_NO_NUMERIC_VALUE
				 && double(numeric_int) == numeric_value
				 && 0 <= numeric_int
				 && numeric_int <= 9);
	  
	    found |= replace;
	    uword_digits.append(replace ? '@' : c);
	  }
	
	  if (found) {
	    uword_digits.append('>');
	  
	    std::string word_digits;
	    uword_digits.toUTF8String(word_digits);
	  
	    pair.second = word_digits;
	  } else
	    pair.second = word;
	}
      }
    
      return pair.second;
    }

  };
};
