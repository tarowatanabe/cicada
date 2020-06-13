//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/hiragana.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>

#include <boost/thread.hpp>

#include <utils/atomicop.hpp>

namespace cicada
{

  namespace stemmer
  {
    static boost::once_flag hiragana_installer_once = BOOST_ONCE_INIT;

    struct HiraganaDetail
    {
      HiraganaDetail()
      {
	boost::call_once(hiragana_installer_once, initialize);
      }
      
      static void initialize()
      {
	// Any-Hiragana, NFKD, remove accents, NFKC
	UErrorCode status = U_ZERO_ERROR;
	UParseError status_parse;
	icu::Transliterator* __trans = icu::Transliterator::createFromRules(icu::UnicodeString::fromUTF8("AnyHiraganaNoAccents"),
									    icu::UnicodeString::fromUTF8(":: Any-Latin; :: NFKD; [[:Z:][:M:][:C:]] > ; :: NFKC; :: Latin-Hiragana;"),
									    UTRANS_FORWARD, status_parse, status);
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("transliterator::create_from_rules(): ") + u_errorName(status));
	
	// register here...
	icu::Transliterator::registerInstance(__trans);
      }
    };
    
    struct HiraganaImpl : public HiraganaDetail
    {
    private:
      icu::Transliterator* trans;
      
    public:
      HiraganaImpl() : HiraganaDetail(), trans(0)
      {
	UErrorCode status = U_ZERO_ERROR;
	trans = icu::Transliterator::createInstance(icu::UnicodeString::fromUTF8("AnyHiraganaNoAccents"),
						    UTRANS_FORWARD,
						    status);
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      }
      ~HiraganaImpl() { std::unique_ptr<icu::Transliterator> tmp(trans); }
      
    public:
      void operator()(icu::UnicodeString& data) { trans->transliterate(data); }
    };
    
    Hiragana::Hiragana() : pimpl(new impl_type()) {}
    Hiragana::~Hiragana() { std::unique_ptr<impl_type> tmp(pimpl); }

    std::string Hiragana::operator()(const utils::piece& word) const
    {
      if (! pimpl)
	throw std::runtime_error("no hiragana normalizer?");

      if (word.empty()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not prefixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      icu::UnicodeString uword = icu::UnicodeString::fromUTF8(icu::StringPiece(word.data(), word.size()));
      
      pimpl->operator()(uword);
      
      if (! uword.isEmpty()) {
	std::string word_hiragana;
	uword.toUTF8String(word_hiragana);
	
	return word_hiragana;
      } else
	return word;
    }

  };
};
