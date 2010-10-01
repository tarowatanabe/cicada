
#include "eval/score.hpp"
#include "eval/per.hpp"
#include "eval/wer.hpp"
#include "eval/ter.hpp"
#include "eval/bleu.hpp"

#include "parameter.hpp"

#include <boost/thread.hpp>

#include "utils/config.hpp"
#include "utils/lexical_cast.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>
#include <unicode/bytestream.h>
#include <unicode/translit.h>

#include <boost/filesystem.hpp>

namespace cicada
{
  namespace eval
  {

    void Scorer::split_non_ascii_characters(const sentence_type& sentence, sentence_type& sentence_split) const
    {
      sentence_split.clear();

      std::string buffer;
	
      sentence_type::const_iterator siter_end = sentence.end();
      for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	  
	UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(*siter));
	  
	StringCharacterIterator iter(uword);
	for (iter.setToStart(); iter.hasNext(); /**/) {
	  const UChar32 c = iter.next32PostInc();
	    
	  if (c < 128)
	    buffer.push_back(c);
	  else {
	    // we will split...
	    if (! buffer.empty())
	      sentence_split.push_back(word_type(buffer.begin(), buffer.end()));
	    buffer.clear();
	      
	    StringByteSink<std::string> __sink(&buffer);
	    UnicodeString(c).toUTF8(__sink);
	      
	    sentence_split.push_back(word_type(buffer.begin(), buffer.end()));
	    buffer.clear();
	  }
	}
	
	if (! buffer.empty())
	  sentence_split.push_back(word_type(buffer.begin(), buffer.end()));
	buffer.clear();
      }
    }
    
    void Scorer::lower_case(const sentence_type& sentence, sentence_type& sentence_lower) const
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

      sentence_lower.clear();
      
      std::string buffer;
      sentence_type::const_iterator siter_end = sentence.end();
      for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(*siter));
	trans.transliterate(uword);
	
	buffer.clear();
	StringByteSink<std::string> __sink(&buffer);
	uword.toUTF8(__sink);

	sentence_lower.push_back(buffer);
      }
    }

    std::string Scorer::lists()
    {
      static const char* desc = "\
bleu:\n\
\torder=<order, default=4> ngram order\n\
\tsplit=[true|false] perform character splitting\n\
per:\n\
\tsplit=[true|false] perform character splitting\n\
wer:\n\
\tsplit=[true|false] perform character splitting\n\
ter:\n\
\tsplit=[true|false] perform character splitting\n\
";

      return desc;
    }
    
    Scorer::scorer_ptr_type Scorer::create(const std::string& parameter)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (param.name() == "bleu" || param.name() == "bleu-linear") {
	int  order = 4;
	bool split = false;
	bool exact = false;

	bool yield_source = false;
	bool yield_target = false;
	
	std::string name;
	path_type   refset_file;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "order") == 0)
	    order = boost::lexical_cast<int>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "exact") == 0)
	    exact = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "split") == 0)
	    split = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "name") == 0)
	    name = piter->second;
	  else if (strcasecmp(piter->first.c_str(), "refset") == 0)
	    refset_file = piter->second;
	  else if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	    const std::string& yield = piter->second;
	    
	    if (strcasecmp(yield.c_str(), "source") == 0)
	      yield_source = true;
	    else if (strcasecmp(yield.c_str(), "target") == 0)
	      yield_target = true;
	    else
	      throw std::runtime_error("unknown parameter: " + parameter);
	    
	  } else
	    std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
	}
	
	return scorer_ptr_type(new BleuScorer(order, split));
      } else if (param.name() == "per") {
	bool split = false;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "split") == 0)
	    split = utils::lexical_cast<bool>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for per: " << piter->first << "=" << piter->second << std::endl;
	}
	
	return scorer_ptr_type(new PERScorer(split));
      } else if (param.name() == "wer") {
	bool split = false;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "split") == 0)
	    split = utils::lexical_cast<bool>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for wer: " << piter->first << "=" << piter->second << std::endl;
	}
	
	return scorer_ptr_type(new WERScorer(split));
      } else if (param.name() == "ter") {
	bool split = false;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "split") == 0)
	    split = utils::lexical_cast<bool>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for ter: " << piter->first << "=" << piter->second << std::endl;
	}
	
	return scorer_ptr_type(new TERScorer(split));
      } else
	throw std::runtime_error("unknown scorer" + param.name());
      
    }
    
  };
};
