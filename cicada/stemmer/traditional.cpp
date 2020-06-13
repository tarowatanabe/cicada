//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/traditional.hpp"

#include <boost/thread.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>

#include "utils/config.hpp"

namespace cicada
{
  namespace stemmer
  {
    Traditional::Traditional() : pimpl(0)
    {
      UErrorCode status = U_ZERO_ERROR;
      std::unique_ptr<icu::Transliterator> trans(icu::Transliterator::createInstance(icu::UnicodeString::fromUTF8("Simplified-Traditional"), UTRANS_FORWARD, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      
      pimpl = trans.release();
    }
    
    Traditional::~Traditional()
    {
      std::unique_ptr<icu::Transliterator> tmp(static_cast<icu::Transliterator*>(pimpl));
    }
    
    std::string Traditional::operator()(const utils::piece& word) const
    {
      if (! pimpl)
	throw std::runtime_error("no traditional?");
      
      if (word.empty()) return word;
      
      const size_type word_size = word.size();
      
      // SGML-like symbols are not traditionaled...
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      icu::UnicodeString uword = icu::UnicodeString::fromUTF8(static_cast<const std::string&>(word));
      
      static_cast<icu::Transliterator*>(pimpl)->transliterate(uword);
      
      std::string word_traditional;
      uword.toUTF8String(word_traditional);
      
      return word_traditional;
    }

  };
};
