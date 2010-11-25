// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TOKENIZER__NONASCII__HPP__
#define __CICADA__TOKENIZER__NONASCII__HPP__ 1

#include <vector>
#include <string>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>
#include <unicode/bytestream.h>

namespace cicada
{
  namespace tokenizer
  {
    class NonAscii : public cicada::Tokenizer
    {
    protected:
      virtual void tokenize(const sentence_type& source, sentence_type& tokenized) const
      {
	std::string buffer;
	
	tokenized.clear();
	sentence_type::const_iterator siter_end = source.end();
	for (sentence_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter) {
	  
	  UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(*siter));
	  
	  StringCharacterIterator iter(uword);
	  for (iter.setToStart(); iter.hasNext(); /**/) {
	    const UChar32 c = iter.next32PostInc();
	    
	    if (c < 128)
	      buffer.push_back(c);
	    else {
	      // we will split...
	      if (! buffer.empty())
		tokenized.push_back(word_type(buffer.begin(), buffer.end()));
	      buffer.clear();
	      
	      StringByteSink<std::string> __sink(&buffer);
	      UnicodeString(c).toUTF8(__sink);
	      
	      tokenized.push_back(word_type(buffer.begin(), buffer.end()));
	      buffer.clear();
	    }
	  }
	  
	  if (! buffer.empty())
	    tokenized.push_back(word_type(buffer.begin(), buffer.end()));
	  buffer.clear();
	}
      }
    };
  };
};

#endif
