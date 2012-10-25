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

    std::string NFKD::operator()(const utils::piece& word) const
    {
      if (word.empty()) return word;
      
      const size_type word_size = word.size();
      
      // SGML-like symbols are not prefixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      icu::UnicodeString uword = icu::UnicodeString::fromUTF8(icu::StringPiece(word.data(), word.size()));
      icu::UnicodeString uword_nfkd;
      
      UErrorCode status = U_ZERO_ERROR;
      static_cast<const icu::Normalizer2*>(handle)->normalize(uword, uword_nfkd, status);
      if (U_FAILURE(status))
	throw std::runtime_error("normalization failed");
      
      std::string word_nfkd;
      uword_nfkd.toUTF8String(word_nfkd);
      
      return word_nfkd;
    }

  };
};
