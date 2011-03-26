//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "grammar_unknown.hpp"

namespace cicada
{
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
    rule_pair_set_type::const_iterator riter_end = __rules.end();
    for (rule_pair_set_type::const_iterator riter = __rules.begin(); riter != riter_end; ++ riter) {
      
      
    }
  }
};
