//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "grammar_unknown.hpp"

#include <utils/lexical_cast.hpp>
#include <utils/space_separator.hpp>

#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>

namespace cicada
{
  void GrammarUnknown::read_character(const std::string& file)
  {
    
    
  }

  void GrammarUnknown::insert(const symbol_type& word)
  {
    id_type node = base_type::next(base_type::root(), word);
    if (node != base_type::root()) return;
    
    // word is oov
    const symbol_type sig = signature->operator()(word);
    node = base_type::next(base_type::root(), sig);
    if (node == base_type::root())
      throw std::runtime_error("invalid signature? " + static_cast<const std::string&>(sig)
			       + " word: " + static_cast<const std::string&>(word));
    
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
      UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
      
      rule_pair_set_type::const_iterator riter_end = __rules.end();
      for (rule_pair_set_type::const_iterator riter = __rules.begin(); riter != riter_end; ++ riter) {
	const rule_ptr_type rule = rule_type::create(rule_type(riter->source->lhs, rule_type::symbol_set_type(1, word)));

	const symbol_type& tag =  rule->lhs;
	
	// compute score by characters...
	double logprob = 0.0;
	StringCharacterIterator iter(uword);
	for (iter.setToStart(); iter.hasNext(); /**/) {
	  const UChar32 ch = iter.next32PostInc();
	  
	  // we will query tag-sig-ch
	  {
	    ngram_set_type::id_type node = ngram.find(ngram.root(), tag);
	    if (node == ngram.root())
	      throw std::runtime_error("cannot find tag in trigram? " + static_cast<const std::string&>(tag));
	    
	    node = ngram.find(node, sig);
	    if (node == ngram.root())
	      throw std::runtime_error("cannot find signature in trigram? " + static_cast<const std::string&>(sig));
	    
	    unigram_set_type::const_iterator titer = ngram[node].find(ch);
	    if (titer != ngram[node].end()) {
	      logprob += titer->second;
	      continue;
	    }
	  }
	  
	  // backoff with tag-sig
	  {
	    backoff_set_type::id_type node = backoff.find(backoff.root(), tag);
	    if (node == backoff.root())
	      throw std::runtime_error("cannot find tag in trigram backoff? " + static_cast<const std::string&>(tag));
	    
	    node = backoff.find(node, sig);
	    if (node == backoff.root())
	      throw std::runtime_error("cannot find signature in trigram backoff? " + static_cast<const std::string&>(sig));
	    
	    logprob += backoff[node];
	  }

	  // we will query sig-ch
	  {
	    ngram_set_type::id_type node = ngram.find(ngram.root(), sig);
	    if (node == ngram.root())
	      throw std::runtime_error("cannot find signature in bigram? " + static_cast<const std::string&>(sig));
	    
	    unigram_set_type::const_iterator titer = ngram[node].find(ch);
	    if (titer != ngram[node].end()) {
	      logprob += titer->second;
	      continue;
	    }
	  }
	  
	  // backoff with sig
	  {
	    backoff_set_type::id_type node = backoff.find(backoff.root(), sig);
	    if (node == backoff.root())
	      throw std::runtime_error("cannot find signature in bigram backoff? " + static_cast<const std::string&>(sig));
	    
	    logprob += backoff[node];
	  }
	  
	  // unigram...
	  {
	    unigram_set_type::const_iterator titer = ngram[node].find(ch);
	    if (titer != ngram[node].end())
	      logprob += titer->second;
	    else
	      logprob += logprob_unk;
	  }
	}
	
	rules_new.push_back(rule_pair_type(rule, riter->target, riter->features, riter->attributes));
	rules_new.back().features[feature_character] = logprob;
      }
    }
    
    rule_pair_set_type::const_iterator riter_end = rules_new.end();
    for (rule_pair_set_type::const_iterator riter = rules_new.begin(); riter != riter_end; ++ riter)
      base_type::insert(*riter);
  }
};
