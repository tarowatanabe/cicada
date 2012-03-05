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
	if (titer->is_terminal() && *titer != vocab_type::BOS && *titer != vocab_type::EOS) {
	  const symbol_type& word = *titer;
	  const size_t pos_scan = (word.id() << 1);
	  const size_t pos_flag = (word.id() << 1) + 1;

	  if (pos_flag >= __non_latin.size())
	    __non_latin.resize(pos_flag + 1, false);
	  
	  if (! __non_latin[pos_scan]) {
	    icu::UnicodeString uword = icu::UnicodeString::fromUTF8(static_cast<const std::string&>(word));
	    
	    bool is_non_latin = false;
	    icu::StringCharacterIterator iter(uword);
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
    
    namespace in {
      struct __attribute_integer : public boost::static_visitor<cicada::AttributeVector::int_type>
      {
	typedef cicada::AttributeVector attribute_set_type;
	
	attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
	attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return 0; }
	attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return 0; }
      };
    };
      
    void InternalNodePenalty::apply_estimate(const edge_type& edge, feature_set_type& features) const
    {
      int internal_node = 0;
      attribute_set_type::const_iterator iter = edge.attributes.find(attr_internal_node);
      if (iter != edge.attributes.end())
	internal_node = boost::apply_visitor(in::__attribute_integer(), iter->second);
      
      if (internal_node)
	features[feature_name()] = - internal_node;
      else
	features.erase(feature_name());
    }
    
  };
};
