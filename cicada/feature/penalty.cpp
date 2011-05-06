#include "penalty.hpp"

#include "utils/lexical_cast.hpp"

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>

namespace cicada
{
  namespace feature
  {
    void ArityPenalty::apply_estimate(const edge_type& edge, feature_set_type& features) const
    {
      size_t count = 0;
      rule_type::symbol_set_type::const_iterator titer_end = edge.rule->rhs.end();
      for (rule_type::symbol_set_type::const_iterator titer = edge.rule->rhs.begin(); titer != titer_end; ++ titer)
	count += titer->is_non_terminal();

      feature_name_set_type& __names = const_cast<feature_name_set_type&>(names);
      
      if (count >= __names.size())
	__names.resize(count + 1);
	
      if (__names[count].empty())
	__names[count] = "arity-penalty:" + utils::lexical_cast<std::string>(count);
	
      features[__names[count]] = -1;
    }

    void NonLatinPenalty::apply_estimate(const edge_type& edge, feature_set_type& features) const
    {
      non_latin_type& __non_latin = const_cast<non_latin_type&>(non_latin);

      int count = 0;
      rule_type::symbol_set_type::const_iterator titer_end = edge.rule->rhs.end();
      for (rule_type::symbol_set_type::const_iterator titer = edge.rule->rhs.begin(); titer != titer_end; ++ titer) 
	if (titer->is_terminal()) {
	  const symbol_type& word = *titer;
	  const size_t pos_scan = (word.id() << 1);
	  const size_t pos_flag = (word.id() << 1) + 1;

	  if (pos_flag >= __non_latin.size())
	    __non_latin.resize(pos_flag + 1, false);
	  
	  if (! __non_latin[pos_scan]) {
	    UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(word));
	    
	    bool is_non_latin = false;
	    StringCharacterIterator iter(uword);
	    for (iter.setToStart(); iter.hasNext(); /**/)
	      is_non_latin |= (iter.next32PostInc() >= 0x80);
	    
	    __non_latin[pos_scan] = true;
	    __non_latin[pos_flag] = is_non_latin;
	  }
	  
	  count += __non_latin[pos_flag];
	}
      
      if (count)
	features[feature_name()] = - count;
      else
	features.erase(feature_name());
    }
  };
};
