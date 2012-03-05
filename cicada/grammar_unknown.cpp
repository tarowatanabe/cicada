//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <string>
#include <vector>

#include "grammar_unknown.hpp"

#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/thread.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>

#include <utils/piece.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/space_separator.hpp>
#include <utils/compress_stream.hpp>

#include "utils/config.hpp"
#include "utils/thread_specific_ptr.hpp"

namespace cicada
{
  

  template <typename Iterator>
  inline
  uint32_t parse_utf8(Iterator first, Iterator last)
  {
    typedef uint32_t uchar;
    
    if (first == last)
      throw std::runtime_error("invalid utf8");
    
    if ((*first & 0x80) == 0)
      return *first;
    else if ((*first & 0xc0) == 0xc0) {
      if (first + 1 >= last)
	throw std::runtime_error("invalid utf8");
      return ((uchar((*first) & 0x1f) << 6)
	      | (uchar((*(first + 1)) & 0x3f)));
    } else if ((*first & 0xe0) == 0xe0) {
      if (first + 2 >= last)
	throw std::runtime_error("invalid utf8");
      return ((uchar((*first) & 0x0f) << 12)
	      | (uchar((*(first + 1)) & 0x3f) << 6)
	      | (uchar((*(first + 2)) & 0x3f)));
    } else if ((*first & 0xf0) == 0xf0) {
      if (first + 3 >= last)
	throw std::runtime_error("invalid utf8");
      return ((uchar((*first) & 0x07) << 18)
	      | (uchar((*(first + 1)) & 0x3f) << 12)
	      | (uchar((*(first + 2)) & 0x3f) << 6)
	      | (uchar((*(first + 3)) & 0x3f)));
    } else if ((*first & 0xf8) == 0xf8) {
            if (first + 4 >= last)
	      throw std::runtime_error("invalid utf8");
      return ((uchar((*first) & 0x03) << 24)
	      | (uchar((*(first + 1)) & 0x3f) << 18)
	      | (uchar((*(first + 2)) & 0x3f) << 12)
	      | (uchar((*(first + 3)) & 0x3f) << 6)
	      | (uchar((*(first + 4)) & 0x3f)));
    } else if ((*first & 0xfc) == 0xfc) {
      if (first + 5 >= last)
	throw std::runtime_error("invalid utf8");
      return ((uchar((*first) & 0x01) << 30)
	      | (uchar((*(first + 1)) & 0x3f) << 24)
	      | (uchar((*(first + 2)) & 0x3f) << 18)
	      | (uchar((*(first + 3)) & 0x3f) << 12)
	      | (uchar((*(first + 4)) & 0x3f) << 6)
	      | (uchar((*(first + 5)) & 0x3f)));
    } else
      throw std::runtime_error("invlaid utf8 sequence");
  }

  template <typename Container>
  inline
  uint32_t parse_utf8(const Container& container)
  {
    return parse_utf8(container.begin(), container.end());
  }
  
  void GrammarUnknown::read_character(const std::string& file)
  {
    typedef std::vector<utils::piece, std::allocator<utils::piece> > tokens_type;
    typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
    
    backoff.clear();
    ngram.clear();
    unigram.clear();
    logprob_unk = 0.0;
    
    utils::compress_istream is(file, 1024 * 1024);
    std::string line;
    tokens_type tokens;
    
    while (std::getline(is, line)) {
      const utils::piece line_piece(line);
      tokenizer_type tokenizer(line_piece);
      
      tokens.clear();
      tokens.insert(tokens.end(), tokenizer.begin(), tokenizer.end());
      
      if (tokens.empty()) continue;

      const double score = utils::lexical_cast<double>(tokens.back());

      if (tokens.front() == "backoff:") {
	switch (tokens.size()) {
	case 3:
	case 4:
	  backoff[backoff.insert(tokens.begin() + 1, tokens.end() - 1)] = score;
	  break;
	default:
	  throw std::runtime_error("invaid backoff? " + line);
	}
      } else {
	switch (tokens.size()) {
	case 5:
	case 4:
	  ngram[ngram.insert(tokens.begin() + 1, tokens.end() - 2)][parse_utf8(*(tokens.end() - 2))] = score;
	  break;
	case 3:
	  if (tokens[1] == "<unk>")
	    logprob_unk = score;
	  else
	    unigram[parse_utf8(tokens[1])] = score;
	  break;
	default:
	  throw std::runtime_error("invaid model? " + line);
	}
      }
    }
  }

  void GrammarUnknown::insert(const symbol_type& word)
  {
    id_type node = base_type::next(base_type::root(), word);
    if (node != base_type::root()) return;
    
    // word is oov
    symbol_type sig = signature->operator()(word);
    node = base_type::next(base_type::root(), sig);
    if (node == base_type::root()) {
      sig =  signature_type::FALLBACK;
      node = base_type::next(base_type::root(), sig);
      
      if (node == base_type::root())
	throw std::runtime_error("invalid signature? " + static_cast<const std::string&>(signature->operator()(word))
				 + " word: " + static_cast<const std::string&>(word));
    }
    
    const rule_pair_set_type& __rules = base_type::rules(node);
    if (__rules.empty())
      throw std::runtime_error("no rules for signature? " + static_cast<const std::string&>(sig)
			       + " word: " + static_cast<const std::string&>(word));
    
    rule_pair_set_type rules_new;
    rules_new.reserve(__rules.size());
    
    if (ngram.empty()) {
      rule_pair_set_type::const_iterator riter_end = __rules.end();
      for (rule_pair_set_type::const_iterator riter = __rules.begin(); riter != riter_end; ++ riter) {
	const rule_ptr_type rule = rule_type::create(rule_type(riter->source->lhs, rule_type::symbol_set_type(1, word)));
	rules_new.push_back(rule_pair_type(rule, riter->target, riter->features, riter->attributes));
      }
    } else {
      icu::UnicodeString uword = icu::UnicodeString::fromUTF8(static_cast<const std::string&>(word));
      
      rule_pair_set_type::const_iterator riter_end = __rules.end();
      for (rule_pair_set_type::const_iterator riter = __rules.begin(); riter != riter_end; ++ riter) {
	const rule_ptr_type rule = rule_type::create(rule_type(riter->source->lhs, rule_type::symbol_set_type(1, word)));

	const symbol_type& tag =  rule->lhs;
	
	// compute score by characters...
	double logprob = 0.0;
	size_t num_char = 0;
	StringCharacterIterator iter(uword);
	for (iter.setToStart(); iter.hasNext(); ++ num_char) {
	  const UChar32 ch = iter.next32PostInc();

	  // we will query sig-tag-ch
	  {
	    ngram_set_type::id_type node = ngram.find(ngram.root(), sig);
	    if (node != ngram.root()) {
	      node = ngram.find(node, tag);
	      
	      if (node != ngram.root()) {
		unigram_set_type::const_iterator titer = ngram[node].find(ch);
		
		if (titer != ngram[node].end()) {
		  logprob += titer->second;
		  continue;
		}
	      }
	    }
	  }

	  // backoff with sig-tag
	  {
	    backoff_set_type::id_type node = backoff.find(backoff.root(), sig);
	    if (node != backoff.root()) {
	      node = backoff.find(node, tag);
	      
	      if (node != backoff.root())
		logprob += backoff[node];
	    }
	  }

	  // we will query tag-ch
	  {
	    ngram_set_type::id_type node = ngram.find(ngram.root(), tag);
	    if (node != ngram.root()) {
	      unigram_set_type::const_iterator titer = ngram[node].find(ch);
	      
	      if (titer != ngram[node].end()) {
		logprob += titer->second;
		continue;
	      }
	    }
	  }

	  // backoff with tag
	  {
	    backoff_set_type::id_type node = backoff.find(backoff.root(), tag);
	    if (node != backoff.root())
	      logprob += backoff[node];
	  }

	  // unigram...
	  {
	    unigram_set_type::const_iterator titer = unigram.find(ch);
	    if (titer != unigram.end())
	      logprob += titer->second;
	    else
	      logprob += logprob_unk;
	  }
	}
	
	rules_new.push_back(rule_pair_type(rule, riter->target, riter->features, riter->attributes));
	rules_new.back().features[feature_character] = logprob / num_char;
      }
    }
    
    rule_pair_set_type::const_iterator riter_end = rules_new.end();
    for (rule_pair_set_type::const_iterator riter = rules_new.begin(); riter != riter_end; ++ riter)
      base_type::insert(*riter);
  }
};
