//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/katakana.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>

#include <boost/thread.hpp>

#include <utils/atomicop.hpp>

namespace cicada
{

  namespace stemmer
  {
    static boost::once_flag katakana_installer_once = BOOST_ONCE_INIT;

    struct KatakanaDetail
    {
      KatakanaDetail()
      {
	boost::call_once(katakana_installer_once, initialize);
      }
      
      static void initialize()
      {
	// Any-Katakana, NFKD, remove accents, NFKC
	UErrorCode status = U_ZERO_ERROR;
	UParseError status_parse;
	icu::Transliterator* __trans = icu::Transliterator::createFromRules(icu::UnicodeString::fromUTF8("AnyKatakanaNoAccents"),
									    icu::UnicodeString::fromUTF8(":: Any-Latin; :: NFKD; [[:Z:][:M:][:C:]] > ; :: NFKC; :: Latin-Katakana;"),
									    UTRANS_FORWARD, status_parse, status);
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("transliterator::create_from_rules(): ") + u_errorName(status));
	
	// register here...
	icu::Transliterator::registerInstance(__trans);
      }
    };
    
    struct KatakanaImpl : public KatakanaDetail
    {
    private:
      icu::Transliterator* trans;
      
    public:
      KatakanaImpl() : KatakanaDetail(), trans(0)
      {
	UErrorCode status = U_ZERO_ERROR;
	trans = icu::Transliterator::createInstance(icu::UnicodeString::fromUTF8("AnyKatakanaNoAccents"),
						    UTRANS_FORWARD,
						    status);
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      }
      ~KatakanaImpl() { std::unique_ptr<icu::Transliterator> tmp(trans); }
      
    public:
      void operator()(icu::UnicodeString& data) { trans->transliterate(data); }
    };
    
    Katakana::Katakana() : pimpl(new impl_type()) {}
    Katakana::~Katakana() { std::unique_ptr<impl_type> tmp(pimpl); }

    std::string Katakana::operator()(const utils::piece& word) const
    {
      if (! pimpl)
	throw std::runtime_error("no katakana normalizer?");

      if (word.empty()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not prefixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      icu::UnicodeString uword = icu::UnicodeString::fromUTF8(icu::StringPiece(word.data(), word.size()));
      
      pimpl->operator()(uword);
      
      if (! uword.isEmpty()) {
	std::string word_katakana;
	uword.toUTF8String(word_katakana);
	
	return word_katakana;
      } else
	return word;
    }

  };
};
