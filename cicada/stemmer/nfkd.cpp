//
//  Copyright(C) 2012-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/nfkd.hpp"

#include <unicode/unistr.h>
#include <unicode/translit.h>

#include <boost/thread.hpp>

#include <utils/atomicop.hpp>

namespace cicada
{

  namespace stemmer
  {
    static boost::once_flag nfkd_installer_once = BOOST_ONCE_INIT;
    
    struct NFKDImpl
    {
      NFKDImpl()
      {
	boost::call_once(nfkd_installer_once, initialize);
      }
      
      static void initialize()
      {
	UErrorCode status = U_ZERO_ERROR;
	UParseError status_parse;
	icu::Transliterator* __trans = icu::Transliterator::createFromRules(icu::UnicodeString::fromUTF8("AnyNFKD"),
									    icu::UnicodeString::fromUTF8("::NFKD; [[:White_Space:]] > ;"),
									    UTRANS_FORWARD, status_parse, status);
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("transliterator::create_from_rules(): ") + u_errorName(status));
	
	// register here...
	icu::Transliterator::registerInstance(__trans);
      }
    };


    NFKD::NFKD()
      : handle(0)
    {
      NFKDImpl __impl;

      UErrorCode status = U_ZERO_ERROR;
      
      std::unique_ptr<icu::Transliterator> trans(icu::Transliterator::createInstance(icu::UnicodeString::fromUTF8("AnyNFKD"), UTRANS_FORWARD, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      
      handle = trans.release();
    }

    NFKD::~NFKD()
    {
      std::unique_ptr<icu::Transliterator> tmp(static_cast<icu::Transliterator*>(handle));
    }

    std::string NFKD::operator()(const utils::piece& word) const
    {
      if (! handle)
	throw std::runtime_error("no NFKD?");
      
      if (word.empty()) return word;
      
      const size_type word_size = word.size();
      
      // SGML-like symbols are not prefixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      icu::UnicodeString uword = icu::UnicodeString::fromUTF8(icu::StringPiece(word.data(), word.size()));
      
      static_cast<icu::Transliterator*>(handle)->transliterate(uword);
      
      std::string word_nfkd;
      uword.toUTF8String(word_nfkd);
      
      return word_nfkd;
    }

  };
};
