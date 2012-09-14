//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <map>

#include <cicada/sentence.hpp>

#include <boost/algorithm/string/erase.hpp>

#include "grammar_format.hpp"

namespace cicada
{
  GrammarFormat::id_type GrammarFormat::next(const id_type& node, const symbol_type& symbol) const
  {
    // check!
    if (symbol.is_non_terminal()) return root();
    
    id_type node_next = base_type::next(node, symbol);
    if (node_next == root()) {
      // the combination of node/symbol has not been visited!
      
      const std::string context = (node == root()
				   ? static_cast<const std::string&>(symbol)
				   : prefix[node] + ' ' + static_cast<const std::string&>(symbol));

      format_type::phrase_set_type phrases;
      if (remove_space)
	format->operator()(boost::algorithm::erase_all_copy(context, " "), phrases);
      else
	format->operator()(context, phrases);
      
      if (! phrases.empty()) {
	typedef cicada::Sentence phrase_type;
	typedef std::map<std::string, feature_set_type, std::less<std::string>, std::allocator<std::pair<const std::string, feature_set_type> > > unique_set_type;

	unique_set_type uniques;
	format_type::phrase_set_type::const_iterator piter_end = phrases.end();
	for (format_type::phrase_set_type::const_iterator piter = phrases.begin(); piter != piter_end; ++ piter)
	  uniques[piter->phrase][static_cast<const std::string&>(feature) + ":" + piter->tag] = 1;
	
	// construct source-side lhs from context
	const phrase_type phrase_source = phrase_type(utils::piece(context));
	const rule_ptr_type rule_source = rule_type::create(rule_type(non_terminal, phrase_source.begin(), phrase_source.end()));
	
	unique_set_type::iterator uiter_end = uniques.end();
	for (unique_set_type::iterator uiter = uniques.begin(); uiter != uiter_end; ++ uiter) {
	  uiter->second[feature] = -1;
	  
	  const phrase_type phrase_target = phrase_type(utils::piece(uiter->first));
	  const rule_ptr_type rule_target = rule_type::create(rule_type(non_terminal, phrase_target.begin(), phrase_target.end()));
	  
	  const_cast<base_type&>(static_cast<const base_type&>(*this)).insert(rule_source, rule_target, uiter->second);
	}
	
	node_next = base_type::next(node, symbol);
      } else
	node_next = const_cast<base_type&>(static_cast<const base_type&>(*this)).insert(node, symbol);
      
      if (node_next >= prefix.size())
	const_cast<prefix_set_type&>(prefix).resize(node_next + 1);
      
      const_cast<prefix_type&>(prefix[node_next]) = context;
    }
    
    return node_next;
  }

};
