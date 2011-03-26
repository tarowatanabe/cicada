//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "grammar_unknown_character.hpp"

#include <boost/filesystem.hpp>

#include "utils/compress_stream.hpp"

namespace cicada
{
  void GrammarUnknownCharacter::read_character(const std::string& path)
  {
    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error("no path? " + path);
    
    
    
  }

  void GrammarUnknownCharacter::insert(const symbol_type& word)
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
    {
      rule_pair_set_type::const_iterator riter_end = __rules.end();
      for (rule_pair_set_type::const_iterator riter = __rules.begin(); riter != riter_end; ++ riter) {
	const rule_ptr_type rule = rule_type::create(rule_type(riter->source->lhs, rule_type::symbol_set_type(1, word)));
	
	// we will add a feature, unknown-character
	
	rules_new.push_back(rule_pair_type(rule, riter->target, riter->features, riter->attributes));
      }
    }
    
    rule_pair_set_type::const_iterator riter_end = rules_new.end();
    for (rule_pair_set_type::const_iterator riter = rules_new.begin(); riter != riter_end; ++ riter)
      base_type::insert(*riter);
  }
};
