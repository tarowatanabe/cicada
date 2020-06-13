//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/lower.hpp"

#include <boost/thread.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>

#include "utils/config.hpp"

namespace cicada
{
  namespace stemmer
  {
    Lower::Lower() : pimpl(0)
    {
      UErrorCode status = U_ZERO_ERROR;
      std::unique_ptr<icu::Transliterator> trans(icu::Transliterator::createInstance(icu::UnicodeString::fromUTF8("Lower"), UTRANS_FORWARD, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      
      pimpl = trans.release();
    }
    
    Lower::~Lower()
    {
      std::unique_ptr<icu::Transliterator> tmp(static_cast<icu::Transliterator*>(pimpl));
    }
    
    std::string Lower::operator()(const utils::piece& word) const
    {
      if (! pimpl)
	throw std::runtime_error("no lower caser?");
      
      if (word.empty()) return word;
      
      const size_type word_size = word.size();
    
      // SGML-like symbols are not lowered...
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      icu::UnicodeString uword = icu::UnicodeString::fromUTF8(icu::StringPiece(word.data(), word.size()));
      
      static_cast<icu::Transliterator*>(pimpl)->transliterate(uword);
      
      std::string word_lower;
      uword.toUTF8String(word_lower);
      
      return word_lower;
    }

  };
};
