#include "stemmer/lower.hpp"

#include <boost/thread.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>
#include <unicode/bytestream.h>

#include "utils/config.hpp"

namespace cicada
{
  namespace stemmer
  {
    Stemmer::symbol_type Lower::operator[](const symbol_type& word) const
    {
#ifdef HAVE_TLS
      static __thread Transliterator* __trans_tls = 0;
      static boost::thread_specific_ptr<Transliterator> __trans;
      if (! __trans_tls) {
	UErrorCode status = U_ZERO_ERROR;
	__trans.reset(Transliterator::createInstance(UnicodeString("Lower", "utf-8"), UTRANS_FORWARD, status));
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
	
	__trans_tls = __trans.get();
      }
      
      Transliterator& trans = *__trans_tls;
#else
      static boost::thread_specific_ptr<Transliterator> __trans;
      if (! __trans.get()) {
	UErrorCode status = U_ZERO_ERROR;
	__trans.reset(Transliterator::createInstance(UnicodeString("Lower", "utf-8"), UTRANS_FORWARD, status));
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      }
      
      Transliterator& trans = *__trans;
#endif

      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
      const size_type word_size = word.size();
    
      // SGML-like symbols are not lowered...
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
      
      if (word.id() >= __cache.size())
	__cache.resize(word.id() + 1, vocab_type::EMPTY);
    
      if (__cache[word.id()] == vocab_type::EMPTY) {
	std::string word_lower;
	UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));

	trans.transliterate(uword);
	
	StringByteSink<std::string> __sink(&word_lower);
	uword.toUTF8(__sink);
	
	__cache[word.id()] = word_lower;
      }
    
      return __cache[word.id()];
    }

  };
};
