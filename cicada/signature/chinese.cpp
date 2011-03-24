// -*- encoding: utf-8 -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "signature/chinese.hpp"

#include <boost/thread.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/regex.h>

#include "utils/config.hpp"

namespace cicada
{
  namespace signature
  {
    
    struct ChineseImpl
    {
      struct Matcher
      {
	Matcher(const char* pattern) : matcher(0) { initialize(pattern); }
	~Matcher() { std::auto_ptr<RegexMatcher> tmp(matcher); }
	
	bool operator()(const UnicodeString& x)
	{
	  matcher->reset(x);
	  
	  UErrorCode status = U_ZERO_ERROR;
	  const bool result = matcher->matches(status);
	  if (U_FAILURE(status))
	    throw std::runtime_error(std::string("RegexMatcher::matches(): ") + u_errorName(status));
	  return result;
	}
	
      private:
	void initialize(const char* pattern)
	{
	  UErrorCode status = U_ZERO_ERROR;
	  matcher = new RegexMatcher(UnicodeString::fromUTF8(pattern), 0, status);
	  if (U_FAILURE(status))
	    throw std::runtime_error(std::string("RegexMatcher: ") + u_errorName(status));
	}
      private:
	RegexMatcher* matcher;
      };

      typedef Matcher matcher_type;
      
      
      ChineseImpl()
	: number_match(".*[0-9０-９一二三四五六七八九十百千万亿零〇○◯].*"),
	  date_match(".*[0-9０-９一二三四五六七八九十百千万亿零〇○◯].*[年月日号]"),
	  ordinal_match("第.*"),
	  proper_name_match(".*[··•․‧∙⋅・].*") {}

      matcher_type number_match;
      matcher_type date_match;
      matcher_type ordinal_match;
      matcher_type proper_name_match;
    };

    Chinese::Chinese() : pimpl(new ChineseImpl()) {}
    Chinese::~Chinese() { std::auto_ptr<ChineseImpl> tmp(static_cast<ChineseImpl*>(pimpl)); }
    
    Signature::symbol_type Chinese::operator[](const symbol_type& word) const
    {
      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
      
      const size_type word_size = word.size();
      
      // SGML-like symbols are used as-is
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
      
      if (word.id() >= __cache.size())
	__cache.resize(word.id() + 1, vocab_type::EMPTY);
      
      if (__cache[word.id()] == vocab_type::EMPTY) {
	typedef ChineseImpl impl_type;

	impl_type& impl = *static_cast<impl_type*>(pimpl);
	
	UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
	
	std::string signature = "UNKNOWN";

	if (impl.date_match(uword))
	  signature += "-DATE";
	else if (impl.number_match(uword)) {
	  signature += "-NUM";
	  if (impl.ordinal_match(uword))
	    signature += "-ORD";
	}
	
	if (impl.proper_name_match(uword))
	  signature += "-PROPER";
	
	__cache[word.id()] = signature;
      }
      
      return __cache[word.id()];
    }
      
  };
};
