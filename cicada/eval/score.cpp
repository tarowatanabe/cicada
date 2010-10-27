
#include "eval/score.hpp"
#include "eval/per.hpp"
#include "eval/wer.hpp"
#include "eval/sb.hpp"
#include "eval/sk.hpp"
#include "eval/ter.hpp"
#include "eval/bleu.hpp"
#include "eval/wlcs.hpp"

#include "stemmer.hpp"
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
      typedef cicada::Stemmer stemmer_type;

#ifdef HAVE_TLS
      static __thread stemmer_type* __stemmer_tls = 0;
      static boost::thread_specific_ptr<stemmer_type> __stemmer;
      if (! __stemmer_tls) {
	__stemmer.reset(&stemmer_type::create("lower"));
	__stemmer_tls = __stemmer.get();
      }
      
      stemmer_type& stemmer = *__stemmer_tls;
#else
      static boost::thread_specific_ptr<stemmer_type> __stemmer;
      if (! __stemmer.get())
	__stemmer.reset(&stemmer_type::create("lower"));
      
      stemmer_type& stemmer = *__stemmer;
#endif

      sentence_lower.clear();
      sentence_type::const_iterator siter_end = sentence.end();
      for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter)
	sentence_lower.push_back(stemmer[*siter]);
    }

    std::string Scorer::lists()
    {
      static const char* desc = "\
bleu:\n\
\torder=<order, default=4> ngram order\n\
\tsplit=[true|false] perform character splitting\n\
\tlower=[true|false] perform lower casing\n\
per: position indenendent error rate\n\
\tsplit=[true|false] perform character splitting\n\
\tlower=[true|false] perform lower casing\n\
wer: word error rate\n\
\tsplit=[true|false] perform character splitting\n\
\tlower=[true|false] perform lower casing\n\
ter: translation error rate\n\
\tsplit=[true|false] perform character splitting\n\
\tlower=[true|false] perform lower casing\n\
sk: string kernel\n\
\tp=order of string kernel (default 4)\n\
\tdecay=decay factor for string kernel (default 0.8)\n\
\tsplit=[true|false] perform character splitting\n\
\tlower=[true|false] perform lower casing\n\
wlcs: weighted longest common subsequence\n\
\talpha=length factor, k^alpha (default 1.0)\n\
\tsplit=[true|false] perform character splitting\n\
\tlower=[true|false] perform lower casing\n\
sb: skip bigram\n\
\twindow=window size (default 4, < 0 for infinity, == 0 for non-skip bigram)\n\
\tsplit=[true|false] perform character splitting\n\
\tlower=[true|false] perform lower casing\n\
";

      return desc;
    }
    
    Scorer::scorer_ptr_type Scorer::create(const std::string& parameter)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);

      scorer_ptr_type scorer;
      
      if (param.name() == "bleu" || param.name() == "bleu-linear") {
	int  order = 4;
	bool split = false;
	bool lower = false;
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
	  else if (strcasecmp(piter->first.c_str(), "lower") == 0)
	    lower = utils::lexical_cast<bool>(piter->second);
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
	
	scorer = scorer_ptr_type(new BleuScorer(order));
	scorer->split = split;
	scorer->lower = lower;
      } else if (param.name() == "per") {
	bool split = false;
	bool lower = false;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "split") == 0)
	    split = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "lower") == 0)
	    lower = utils::lexical_cast<bool>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for per: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new PERScorer());
	scorer->split = split;
	scorer->lower = lower;
      } else if (param.name() == "wer") {
	bool split = false;
	bool lower = false;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "split") == 0)
	    split = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "lower") == 0)
	    lower = utils::lexical_cast<bool>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for wer: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new WERScorer());
	scorer->split = split;
	scorer->lower = lower;
      } else if (param.name() == "ter") {
	bool split = false;
	bool lower = false;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "split") == 0)
	    split = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "lower") == 0)
	    lower = utils::lexical_cast<bool>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for ter: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new TERScorer());
	scorer->split = split;
	scorer->lower = lower;

      } else if (param.name() == "sk") {
	int p = 4;
	double decay = 0.8;
	
	bool split = false;
	bool lower = false;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "split") == 0)
	    split = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "lower") == 0)
	    lower = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "p") == 0)
	    p = boost::lexical_cast<int>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "decay") == 0)
	    decay = boost::lexical_cast<double>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for sk: " << piter->first << "=" << piter->second << std::endl;
	}

	scorer = scorer_ptr_type(new SKScorer(p, decay));
	scorer->split = split;
	scorer->lower = lower;
      } else if (param.name() == "wlcs") {
	double alpha = 1.0;
	
	bool split = false;
	bool lower = false;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "split") == 0)
	    split = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "lower") == 0)
	    lower = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "alpha") == 0)
	    alpha = boost::lexical_cast<double>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for wlcs: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new WLCSScorer(alpha));
	scorer->split = split;
	scorer->lower = lower;
	
      } else if (param.name() == "sb") {
	int window = 4;
	
	bool split = false;
	bool lower = false;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "split") == 0)
	    split = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "lower") == 0)
	    lower = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "window") == 0)
	    window = boost::lexical_cast<int>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for sb: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new SBScorer(window));
	scorer->split = split;
	scorer->lower = lower;
	
      } else
	throw std::runtime_error("unknown scorer" + param.name());
      
      return scorer;
    }
    
  };
};
