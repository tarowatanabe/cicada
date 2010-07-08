
#include <utility>
#include <memory>

#include "cicada/feature/ngram_tree.hpp"
#include "cicada/parameter.hpp"

#include "utils/indexed_set.hpp"
#include "utils/compact_trie.hpp"

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace feature
  {
    
    class NGramTreeImpl
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
      
      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;
      
      typedef uint32_t id_type;
      
      struct string_hash_type : public utils::hashmurmur<size_t>
      {
	typedef utils::hashmurmur<size_t> hasher_type;
	
	size_t operator()(const std::string& x) const
	{
	  return hasher_type::operator()(x.begin(), x.end(), 0);
	}
      };
      
      typedef utils::indexed_set<std::string, string_hash_type, std::equal_to<std::string>, std::allocator<std::string> > tree_map_type;
      
      virtual ~NGramTreeImpl() {}
      
      virtual void ngram_tree_score(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features) const = 0;
      virtual void ngram_tree_final_score(const state_ptr_type& state,
					  feature_set_type& features) const = 0;
      
      void clear()
      {
	tree_map.clear();
      }
      
      tree_map_type  tree_map;
      
      phrase_span_set_type phrase_spans_impl;

      bool forced_feature;
    };
    
    template <typename Extract>
    class __NGramTreeImpl : public NGramTreeImpl, public Extract
    {
      
      virtual void ngram_tree_score(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features) const
      {
	// this feature function is complicated in that we know nothing about the source-side...
	
	const std::string& epsilon = static_cast<const std::string&>(vocab_type::EPSILON);
	
	const rule_type::symbol_set_type& phrase = extract_phrase(edge);
	
	if (states.empty()) {
	  // we do not add feature here, since we know nothing abount surrounding context...
	  std::string prefix = epsilon;
	  std::string suffix = epsilon;
	  
	  compute_bound(features, edge, phrase.begin(), phrase.end(), prefix, suffix);
	  
	  id_type* context = reinterpret_cast<id_type*>(state);
	  context[0] = tree_id(compose_path(edge.rule->lhs, prefix));
	  context[1] = tree_id(compose_path(edge.rule->lhs, suffix));
	} else {
	  phrase_span_set_type& phrase_spans = const_cast<phrase_span_set_type&>(phrase_spans_impl);
	  
	  phrase_spans.clear();
	  phrase.terminals(std::back_inserter(phrase_spans));

	  if (phrase_spans.size() != states.size() + 1)
	    throw std::runtime_error("# of states does not match...");
	  
	  std::string prefix;
	  std::string suffix;
	  
	  compute_bound(features, edge, phrase_spans.front().first, phrase_spans.front().second, prefix, suffix);
	  
	  phrase_span_set_type::const_iterator siter_begin = phrase_spans.begin();
	  phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	  for (phrase_span_set_type::const_iterator siter = siter_begin + 1; siter != siter_end; ++ siter) {
	    const phrase_span_type& span = *siter;
	    
	    // incase, we are working with non-synchronous parsing!
	    int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      antecedent_index = siter - (siter_begin + 1);
	    
	    const id_type* antecedent_context = reinterpret_cast<const id_type*>(states[antecedent_index]);
	    
	    const std::string prefix_antecedent = tree_map[antecedent_context[0]];
	    const std::string suffix_antecedent = tree_map[antecedent_context[1]];
	    
	    std::string prefix_next;
	    std::string suffix_next;
	    compute_bound(features, edge, span.first, span.second, prefix_next, suffix_next);
	    
	    if (! suffix.empty())
	      apply_feature(features, edge.rule->lhs, suffix, prefix_antecedent);
	    
	    if (prefix.empty())
	      prefix = prefix_antecedent;
	    
	    if (! prefix_next.empty()) {
	      apply_feature(features, edge.rule->lhs, suffix_antecedent, prefix_next);
	      suffix = suffix_next;
	    } else
	      suffix = suffix_antecedent;
	  }
	  
	  // construct state...
	  id_type* context = reinterpret_cast<id_type*>(state);
	  context[0] = tree_id(compose_path(edge.rule->lhs, prefix.empty() ? epsilon : prefix));
	  context[1] = tree_id(compose_path(edge.rule->lhs, suffix.empty() ? epsilon : suffix));
	}
      }

      id_type tree_id(const std::string& node) const
      {
	tree_map_type::iterator iter = const_cast<tree_map_type&>(tree_map).insert(node).first;
	return iter - tree_map.begin();
      }

      void apply_feature(feature_set_type& features, const std::string& node, const std::string& prev, const std::string& next) const
      {
	const std::string name = feature_name(node, prev, next);
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;
      }

      template <typename Iterator>
      void compute_bound(feature_set_type& features, const edge_type& edge, Iterator first, Iterator last, std::string& prefix, std::string& suffix) const
      {
	bool is_front = true;
	for (/**/; first != last; ++ first) 
	  if (*first != vocab_type::EPSILON) {
	    if (is_front)
	      prefix = *first;
	    
	    // we will count only the boundary...
	    //if (suffix != epsilon)
	    //  apply_feature(features, edge.rule->lhs, suffix, *first);
	    
	    suffix = *first;
	    is_front = false;
	  }
      }
      
      
      const std::string compose_path(const std::string& node, const std::string& antecedent) const
      {
	return node + '(' + antecedent + ')';
      }

      const std::string compose_tree(const std::string& node, const std::string& prev, const std::string& next) const
      {
	return node + "(" + prev + ")(" + next + ')';
      }

      const std::string feature_name(const std::string& node, const std::string& prev, const std::string& next) const
      {
	return Extract::feature_prefix +  compose_tree(node, prev, next);
      }

      
      virtual void ngram_tree_final_score(const state_ptr_type& state,
					  feature_set_type& features) const
      {
	const id_type* antecedent_context = reinterpret_cast<const id_type*>(state);
	
	const std::string prefix_antecedent = tree_map[antecedent_context[0]];
	const std::string suffix_antecedent = tree_map[antecedent_context[1]];
	
	apply_feature(features, vocab_type::GOAL, vocab_type::BOS, prefix_antecedent);
	apply_feature(features, vocab_type::GOAL, vocab_type::EOS, suffix_antecedent);
      }
	  

      template <typename Edge>
      const rule_type::symbol_set_type& extract_phrase(const Edge& x) const
      {
	static const rule_type::symbol_set_type __tmptmp;

	return Extract::operator()(x, __tmptmp);
      }
    };
    

    struct __ngram_tree_extract_source
    {
      __ngram_tree_extract_source()
	: feature_prefix("ngram-tree-source:") {}
      
      template <typename Edge, typename Phrase>
      const Phrase& operator()(const Edge& x, const Phrase& phrase) const
      {
	return x.rule->source;
      }
      
      const std::string feature_prefix;
    };

    struct __ngram_tree_extract_target
    {
      __ngram_tree_extract_target()
	: feature_prefix("ngram-tree-target:") {}
      
      template <typename Edge, typename Phrase>
      const Phrase& operator()(const Edge& x, const Phrase& phrase) const
      {
	return x.rule->target;
      }
      
      const std::string feature_prefix;
    };
    
    
    NGramTree::NGramTree(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "ngram-tree")
	throw std::runtime_error("is this really ngram tree feature function? " + parameter);

      bool source = false;
      bool target = false;
      if (param.find("yield") != param.end()) {
	const std::string& yield = param.find("yield")->second;

	if (strcasecmp(yield.c_str(), "source") == 0)
	  source = true;
	else if (strcasecmp(yield.c_str(), "target") == 0)
	  target = true;
	else
	  throw std::runtime_error("unknown parameter: " + parameter);
      }
      
      if (source && target)
	throw std::runtime_error("both source and target?");
      if (! source && ! target)
	throw std::runtime_error("what side are you going to use?");
      
      std::auto_ptr<impl_type> ngram_tree_impl(source
					       ? dynamic_cast<impl_type*>(new __NGramTreeImpl<__ngram_tree_extract_source>())
					       : dynamic_cast<impl_type*>(new __NGramTreeImpl<__ngram_tree_extract_target>()));

      
      // non-terminal + two neighbouring symbols + span-size
      base_type::__state_size = sizeof(impl_type::id_type) * 2;
      base_type::__feature_name = std::string("ngram-tree-") + (source ? "source" : "target");
      base_type::__sparse_feature = true;
      
      pimpl = ngram_tree_impl.release();
    }
    
    NGramTree::~NGramTree() { std::auto_ptr<impl_type> tmp(pimpl); }

    template <typename FeaturePrefix, typename Feature>
    inline
    bool equal_prefix(const FeaturePrefix& prefix, const Feature& x)
    {
      return x.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), x.begin());
    }
    
    void NGramTree::operator()(state_ptr_type& state,
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

      pimpl->ngram_tree_score(state, states, edge, features);
    }
    
    void NGramTree::operator()(const state_ptr_type& state,
			       feature_set_type& features,
			       feature_set_type& estimates) const
    {
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();

      pimpl->ngram_tree_final_score(state, features);
    }

    void NGramTree::initialize()
    {
      pimpl->clear();
    }
  };
};
