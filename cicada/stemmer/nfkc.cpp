//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/nfkc.hpp"

#include <unicode/unistr.h>
#include <unicode/normalizer2.h>

namespace cicada
{

  namespace stemmer
  {
    NFKC::NFKC()
      : handle(0)
    {
      UErrorCode status = U_ZERO_ERROR;
      handle = icu::Normalizer2::getInstance(NULL, "nfkc", UNORM2_COMPOSE, status);
    }

    std::string NFKC::operator()(const utils::piece& word) const
    {
      if (word.empty()) return word;
      
      const size_type word_size = word.size();
      
      // SGML-like symbols are not prefixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      icu::UnicodeString uword = icu::UnicodeString::fromUTF8(icu::StringPiece(word.data(), word.size()));
      icu::UnicodeString uword_nfkc;
      
      UErrorCode status = U_ZERO_ERROR;
      static_cast<const icu::Normalizer2*>(handle)->normalize(uword, uword_nfkc, status);
      if (U_FAILURE(status))
	throw std::runtime_error("normalization failed");
      
      std::string word_nfkc;
      uword_nfkc.toUTF8String(word_nfkc);
      
      return word_nfkc;
    }

  };
};
