//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/nfkd.hpp"

#include <unicode/unistr.h>
#include <unicode/normalizer2.h>

namespace cicada
{

  namespace stemmer
  {
    NFKD::NFKD()
      : handle(0)
    {
      UErrorCode status = U_ZERO_ERROR;
      handle = icu::Normalizer2::getInstance(NULL, "nfkd", UNORM2_COMPOSE, status);
    }


    Stemmer::symbol_type NFKD::operator[](const symbol_type& word) const
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
	icu::UnicodeString uword_nfkd;
	
	UErrorCode status = U_ZERO_ERROR;
	static_cast<const icu::Normalizer2*>(handle)->normalize(uword, uword_nfkd, status);
	if (U_FAILURE(status))
	  throw std::runtime_error("normalization failed");
	
	std::string word_nfkd;
	uword_nfkd.toUTF8String(word_nfkd);
	
	pair.first  = word;
	pair.second = word_nfkd;
      }
      return pair.second;
    }

  };
};
