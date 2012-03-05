//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer/latin.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>

#include <boost/thread.hpp>

#include <utils/atomicop.hpp>

namespace cicada
{

  namespace stemmer
  {
    static boost::once_flag latin_installer_once = BOOST_ONCE_INIT;

    struct LatinDetail
    {
      LatinDetail()
      {
	boost::call_once(latin_installer_once, initialize);
      }
      
      static void initialize()
      {
	// Any-Latin, NFKD, remove accents, NFKC
	UErrorCode status = U_ZERO_ERROR;
	UParseError status_parse;
	icu::Transliterator* __trans = icu::Transliterator::createFromRules(icu::UnicodeString::fromUTF8("AnyLatinNoAccents"),
									    icu::UnicodeString::fromUTF8(":: Any-Latin; :: NFKD; [[:Z:][:M:][:C:]] > ; :: NFKC;"),
									    UTRANS_FORWARD, status_parse, status);
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("transliterator::create_from_rules(): ") + u_errorName(status));
	
	// register here...
	icu::Transliterator::registerInstance(__trans);
      }
    };
    
    struct LatinImpl : public LatinDetail
    {
    private:
      icu::Transliterator* trans;
      
    public:
      LatinImpl() : LatinDetail(), trans(0)
      {
	UErrorCode status = U_ZERO_ERROR;
	trans = icu::Transliterator::createInstance(icu::UnicodeString::fromUTF8("AnyLatinNoAccents"),
						    UTRANS_FORWARD,
						    status);
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      }
      ~LatinImpl() { std::auto_ptr<icu::Transliterator> tmp(trans); }
      
    public:
      void operator()(icu::UnicodeString& data) { trans->transliterate(data); }
    };
    
    Latin::Latin() : pimpl(new impl_type()) {}
    Latin::~Latin() { std::auto_ptr<impl_type> tmp(pimpl); }

    Stemmer::symbol_type Latin::operator[](const symbol_type& word) const
    {
      if (! pimpl)
	throw std::runtime_error("no latin normalizer?");

      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not prefixed
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
    
      symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
    
      if (word.id() >= __cache.size())
	__cache.resize(word.id() + 1, vocab_type::EMPTY);
      
      if (__cache[word.id()] == vocab_type::EMPTY) {
	
	icu::UnicodeString uword = icu::UnicodeString::fromUTF8(static_cast<const std::string&>(word));
	
	pimpl->operator()(uword);
      
	if (! uword.isEmpty()) {
	  std::string word_latin;
	  uword.toUTF8String(word_latin);
	
	  __cache[word.id()] = word_latin;
	} else
	  __cache[word.id()] = word;
      }
    
      return __cache[word.id()];
    }

  };
};
