
#include <utility>
#include <memory>

#include "cicada/feature/parent.hpp"
#include "cicada/parameter.hpp"

#include "utils/indexed_set.hpp"
#include "utils/compact_trie.hpp"
#include "utils/lexical_cast.hpp"

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace feature
  {
    
    class ParentImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef rule_type::symbol_set_type phrase_type;
      
      
      struct string_hash : public utils::hashmurmur<size_t>
      {
	typedef utils::hashmurmur<size_t> hasher_type;
	
	size_t operator()(const std::string& x) const
	{
	  return hasher_type::operator()(x.begin(), x.end(), 0);
	}
      };
      
      typedef utils::indexed_set<std::string, string_hash, std::equal_to<std::string>, std::allocator<std::string> > string_map_type;
      
      typedef string_map_type::index_type id_type;
      
      ParentImpl()
	: forced_feature(false) {}

      virtual ~ParentImpl() {}
      
      virtual void parent_score(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features) const = 0;
      virtual void parent_final_score(const state_ptr_type& state,
				      const edge_type& edge,
				      feature_set_type& features) const = 0;
      
      void clear()
      {
	string_map.clear();
      }
      
      string_map_type string_map;
      
      bool forced_feature;
    };
    
    template <typename Extract>
    class __ParentImpl : public ParentImpl, public Extract
    {
      
      virtual void parent_score(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features) const
      {
	// this feature function is complicated in that we know nothing about the source-side...
	
	const rule_type::symbol_set_type& phrase = extract_phrase(edge);
	
	if (states.empty()) {
	  std::string rule_string = edge.rule->lhs;
	  if (! phrase.empty()) {
	    rule_string += '(';
	    phrase_type::const_iterator piter_end = phrase.end();
	    for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end - 1; ++ piter)
	      rule_string += static_cast<const std::string&>(*piter) + '|';
	    rule_string += static_cast<const std::string&>(*(piter_end - 1)) + ')';
	  }
	  
	  apply_feature(features, rule_string);
	  
	  string_map_type::iterator iter = const_cast<string_map_type&>(string_map).insert(rule_string).first;
	  
	  *reinterpret_cast<id_type*>(state) = iter - string_map.begin();
	} else {
	  std::string rule_string;
	  
	  int pos_non_terminal = 0;
	  phrase_type::const_iterator piter_end = phrase.end();
	  for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter)
	    if (piter->is_non_terminal()) {
	      int antecedent_index = piter->non_terminal_index() - 1;
	      if (antecedent_index < 0)
		antecedent_index = pos_non_terminal;
	      
	       const id_type* antecedent = reinterpret_cast<const id_type*>(states[antecedent_index]);
	      
	      apply_feature(features, edge.rule->lhs, string_map[*antecedent]);
	      
	      if (rule_string.empty())
		rule_string += static_cast<const std::string&>(piter->non_terminal());
	      else
		rule_string += '|' + static_cast<const std::string&>(piter->non_terminal());
	      
	      ++ pos_non_terminal;
	    } else {
	      if (rule_string.empty())
		rule_string += static_cast<const std::string&>(*piter);
	      else
		rule_string += '|' + static_cast<const std::string&>(*piter);
	    }
	  
	  rule_string = static_cast<const std::string&>(edge.rule->lhs) + '(' + rule_string + ')';
	  
	  apply_feature(features, rule_string);
	  
	  string_map_type::iterator iter = const_cast<string_map_type&>(string_map).insert(rule_string).first;
	  
	  *reinterpret_cast<id_type*>(state) = iter - string_map.begin();
	}
      }
      
      virtual void parent_final_score(const state_ptr_type& state,
				      const edge_type& edge,
				      feature_set_type& features) const
      {
	const id_type* id = reinterpret_cast<const id_type*>(state);
	
	apply_feature(features, vocab_type::BOS, string_map[*id]);
      }
      
      void apply_feature(feature_set_type& features,
			 const std::string& rule) const
      {
	const std::string name = Extract::feature_prefix + rule;
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;
      }
      
      void apply_feature(feature_set_type& features,
			 const std::string& parent,
			 const std::string& rule) const
      {
	const std::string name = Extract::feature_prefix + parent + '(' + rule + ')';
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;
      }
      

      template <typename Edge>
      const rule_type::symbol_set_type& extract_phrase(const Edge& x) const
      {
	static const rule_type::symbol_set_type __tmptmp;

	return Extract::operator()(x, __tmptmp);
      }
    };
    

    struct __parent_extract_source
    {
      __parent_extract_source()
	: feature_prefix("parent-source:") {}
      
      template <typename Edge, typename Phrase>
      const Phrase& operator()(const Edge& x, const Phrase& phrase) const
      {
	return x.rule->source;
      }
      
      const std::string feature_prefix;
    };

    struct __parent_extract_target
    {
      __parent_extract_target()
	: feature_prefix("parent-target:") {}
      
      template <typename Edge, typename Phrase>
      const Phrase& operator()(const Edge& x, const Phrase& phrase) const
      {
	return x.rule->target;
      }
      
      const std::string feature_prefix;
    };

    
    Parent::Parent(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "parent")
	throw std::runtime_error("is this really parent feature function? " + parameter);

      bool source = false;
      bool target = false;
      
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	  const std::string& yield = piter->second;
	  
	  if (strcasecmp(yield.c_str(), "source") == 0)
	    source = true;
	  else if (strcasecmp(yield.c_str(), "target") == 0)
	    target = true;
	  else
	    throw std::runtime_error("unknown parameter: " + parameter);
	} else
	  std::cerr << "WARNING: unsupported parameter for parent: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (source && target)
	throw std::runtime_error("both source and target?");
      if (! source && ! target)
	throw std::runtime_error("what side are you going to use?");

      std::auto_ptr<impl_type> parent_impl(source
					   ? dynamic_cast<impl_type*>(new __ParentImpl<__parent_extract_source>())
					   : dynamic_cast<impl_type*>(new __ParentImpl<__parent_extract_target>()));
      
      
      // parent conext + terminal-boundary + span-size
      base_type::__state_size = sizeof(impl_type::id_type) + sizeof(symbol_type) * 2 + sizeof(int);
      base_type::__feature_name = std::string("parent-") + (source ? "source" : "target");
      base_type::__sparse_feature = true;
      
      pimpl = parent_impl.release();
    }
    
    Parent::~Parent() { std::auto_ptr<impl_type> tmp(pimpl); }

    
    Parent::Parent(const Parent& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(0)
    {
      typedef __ParentImpl<__parent_extract_source> parent_source_type;
      typedef __ParentImpl<__parent_extract_target> parent_target_type;
      
      if (dynamic_cast<const parent_source_type*>(x.pimpl))
	pimpl = new parent_source_type();
      else
	pimpl = new parent_target_type();
    }
    
    Parent& Parent::operator=(const Parent& x)
    {
      typedef __ParentImpl<__parent_extract_source> parent_source_type;
      typedef __ParentImpl<__parent_extract_target> parent_target_type;

      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      std::auto_ptr<impl_type> tmp(pimpl);
      
      if (dynamic_cast<const parent_source_type*>(x.pimpl))
	pimpl = new parent_source_type();
      else
	pimpl = new parent_target_type();
      
      return *this;
    }


    template <typename FeaturePrefix, typename Feature>
    inline
    bool equal_prefix(const FeaturePrefix& prefix, const Feature& x)
    {
      return x.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), x.begin());
    }
    
    void Parent::operator()(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       feature_set_type& estimates) const
    {
      const std::string& __feature_prefix = base_type::feature_name();
      for (feature_set_type::iterator fiter = features.begin(); fiter != features.end(); /**/)
	if (equal_prefix(__feature_prefix, fiter->first))
	  features.erase(fiter ++);
	else
	  ++ fiter;
      
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();

      pimpl->parent_score(state, states, edge, features);
    }
    
    void Parent::operator()(const state_ptr_type& state,
			    const edge_type& edge,
			    feature_set_type& features,
			    feature_set_type& estimates) const
    {
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->parent_final_score(state, edge, features);
    }

    void Parent::initialize()
    {
      pimpl->clear();
    }
  };
};
