// -*- mode: c++ -*-

#ifndef __CICADA__TOKENIZE__NONASCII__HPP__
#define __CICADA__TOKENIZE__NONASCII__HPP__ 1

#include <vector>
#include <string>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>
#include <unicode/bytestream.h>

namespace cicada
{
  namespace tokenize
  {
    template <typename Sent, typename Tokenized>
    void nonascii(const Sent& sentence, Tokenized& __tokenized)
    {
      typedef Sent sentence_type;
      typedef typename sentence_type::value_type word_type;
      
      std::string buffer;
      
      std::vector<word_type, std::allocator<word_type> > tokenized;
      
      typename sentence_type::const_iterator siter_end = sentence.end();
      for (typename sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	  
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

      __tokenized = Tokenized(tokenized.begin(), tokenized.end());
    }
    
  };
};

#endif
