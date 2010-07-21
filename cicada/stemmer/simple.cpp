
#include "stemmer_simple.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>
#include <unicode/bytestream.h>
#include <unicode/translit.h>
#include <unicode/regex.h>

namespace cicada
{
  
  Stemmer::symbol_type StemmerPrefix::operator[](const symbol_type& word) const
  {
    if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
    const size_type word_size = word.size();
    
    // SGML-like symbols are not prefixed
    if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
      return word;

    symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
    
    if (word.id() >= __cache.size())
      __cache.resize(word.id() + 1, vocab_type::EMPTY);

    if (__cache[word.id()] == vocab_type::EMPTY) {
      UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
    
      const size_t index = uword.moveIndex32(0, int(size));
      
      UnicodeString uword_prefix;
      uword.extractBetween(0, index, uword_prefix);
      
      if (uword_prefix.length() < uword.length()) {
	uword_prefix.append('+');
	
	std::string word_prefix;
	StringByteSink<std::string> __sink(&word_prefix);
	uword_prefix.toUTF8(__sink);
	
	__cache[word.id()] = word_prefix;
      } else
	__cache[word.id()] = word;
    }
    
    return __cache[word.id()];
  }
  
  Stemmer::symbol_type StemmerSuffix::operator[](const symbol_type& word) const
  {
    if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
    const size_type word_size = word.size();
    
    // SGML-like symbols are not suffixed
    if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
      return word;

    symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
    
    if (word.id() >= __cache.size())
      __cache.resize(word.id() + 1, vocab_type::EMPTY);
    
    if (__cache[word.id()] == vocab_type::EMPTY) {
      UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
      
      const size_t index = uword.moveIndex32(uword.length(), - int(size));
      
      UnicodeString uword_suffix;
      uword.extractBetween(index, uword.length(), uword_suffix);
      
      if (uword_suffix.length () < uword.length()) {
	uword_suffix.insert(0, '+');
	
	std::string word_suffix;
	StringByteSink<std::string> __sink(&word_suffix);
	uword_suffix.toUTF8(__sink);
	
	__cache[word.id()] = word_suffix;
      } else
	__cache[word.id()] = word;
    }
    
    return __cache[word.id()];
  }

  Stemmer::symbol_type StemmerDigit::operator[](const symbol_type& word) const
  {
    if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
    const size_type word_size = word.size();
    
    // SGML-like symbols are not prefixed
    if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
      return word;
    
    symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
    
    if (word.id() >= __cache.size())
      __cache.resize(word.id() + 1, vocab_type::EMPTY);
    
    if (__cache[word.id()] == vocab_type::EMPTY) {
      
      UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
      
      if (u_getIntPropertyValue(uword.char32At(0), UCHAR_NUMERIC_TYPE) == U_NT_NONE
	  && u_getIntPropertyValue(uword.char32At(uword.length() - 1), UCHAR_NUMERIC_TYPE) == U_NT_NONE)
	__cache[word.id()] = word;
      else {
	bool found = false;
	UnicodeString uword_digits("<digit-");
	StringCharacterIterator iter(uword);
	for (iter.setToStart(); iter.hasNext(); /**/) {
	  UChar32 c = iter.next32PostInc();
	  
	  //const int32_t numeric_type = u_getIntPropertyValue(c, UCHAR_NUMERIC_TYPE);
	  const double numeric_value = u_getNumericValue(c);
	  const int32_t numeric_int = int(numeric_value);
	  
	  const bool replace =(numeric_value != U_NO_NUMERIC_VALUE
			       && double(numeric_int) == numeric_value
			       && 0 <= numeric_int
			       && numeric_int <= 9);
	  
	  found |= replace;
	  uword_digits.append(replace ? '@' : c);
	}
	
	if (found) {
	  uword_digits.append('>');
	  
	  std::string word_digits;
	  StringByteSink<std::string> __sink(&word_digits);
	  uword_digits.toUTF8(__sink);
	  
	  __cache[word.id()] = word_digits;
	} else
	  __cache[word.id()] = word;
      }
    }
    
    return __cache[word.id()];
  }

  class AnyLatin
  {
  private:
    Transliterator* trans;
    
  public:
    AnyLatin() : trans(0) { open(); }
    ~AnyLatin() { close(); }
    
  private:
    AnyLatin(const AnyLatin& x) { }
    AnyLatin& operator=(const AnyLatin& x) { return *this; }
    
  public:
    void operator()(UnicodeString& data) { trans->transliterate(data); }
    
  private:
    void open()
    {
      close();
      
      __initialize();
      
      UErrorCode status = U_ZERO_ERROR;
      trans = Transliterator::createInstance(UnicodeString::fromUTF8("AnyLatinNoAccents"),
					     UTRANS_FORWARD,
					     status);
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
      
    }
    
    void close()
    {
      if (trans) 
	delete trans;
      trans = 0;
    }
    
  private:
    void __initialize()
    {
      static bool __initialized = false;
      
      if (__initialized) return;
      
      // Any-Latin, NFKD, remove accents, NFKC
      UErrorCode status = U_ZERO_ERROR;
      UParseError status_parse;
      Transliterator* __trans = Transliterator::createFromRules(UnicodeString::fromUTF8("AnyLatinNoAccents"),
								UnicodeString::fromUTF8(":: Any-Latin; :: NFKD; [[:Z:][:M:][:C:]] > ; :: NFKC;"),
								UTRANS_FORWARD, status_parse, status);
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("transliterator::create_from_rules(): ") + u_errorName(status));
      
      // register here...
      Transliterator::registerInstance(__trans);
      
      __initialized = true;
    }
  };

  Stemmer::symbol_type StemmerLatin::operator[](const symbol_type& word) const
  {
    if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
    const size_type word_size = word.size();
    
    // SGML-like symbols are not prefixed
    if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
      return word;
    
    symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
    
    if (word.id() >= __cache.size())
      __cache.resize(word.id() + 1, vocab_type::EMPTY);
    
    if (__cache[word.id()] == vocab_type::EMPTY) {
#ifdef HAVE_TLS
      static __thread AnyLatin* __any_latin_tls = 0;
      static boost::thread_specific_ptr<AnyLatin> __any_latin_specific;
      
      if (! __any_latin_tls) {
	__any_latin_specific.reset(new AnyLatin());
	__any_latin_tls = __any_latin_specific.get();
      }
      
      AnyLatin& any_latin = *__any_latin_tls;
#else
      static boost::thread_specific_ptr<AnyLatin> __any_latin_specific;
      
      if (! __any_latin_specific.get())
	__any_latin_specific.reset(new AnyLatin());
      
      AnyLatin& any_latin = *__any_latin_specific;
#endif
      
      UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
      
      any_latin(uword);
      
      if (! uword.isEmpty()) {
	std::string word_latin;
	StringByteSink<std::string> __sink(&word_latin);
	uword.toUTF8(__sink);
	
	__cache[word.id()] = word_latin;
      } else
	__cache[word.id()] = word;
    }
    
    return __cache[word.id()];
  }

};
