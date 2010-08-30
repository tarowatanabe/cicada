
#include <utility>
#include <memory>

#include "cicada/feature/parent.hpp"
#include "cicada/parameter.hpp"

#include "utils/indexed_set.hpp"
#include "utils/lexical_cast.hpp"

#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"


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
      typedef cicada::Cluster  cluster_type;
      typedef cicada::Stemmer  stemmer_type;
      
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
	: cluster(0), stemmer_prefix(0), stemmer_suffix(0), stemmer_digits(0), forced_feature(false) {}

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

      cluster_type* cluster;
      stemmer_type* stemmer_prefix;
      stemmer_type* stemmer_suffix;
      stemmer_type* stemmer_digits;
      
      string_map_type string_map;
      
      bool forced_feature;
    };

    struct __extractor_none
    {
      template <typename Word>
      std::string operator()(const Word& word) const
      {
	return word;
      }
    };

    template <typename Extract>
    struct __extractor
    {
      __extractor(const Extract& __extract) : extract(__extract) {}
      
      const Extract& extract;
      
      template <typename Word>
      std::string operator()(const Word& word) const
      {
	return extract[word];
      }
    };
    
    template <typename Iterator, typename Extractor>
    std::string extract_phrase_rule(const std::string& lhs,
				    Iterator first, Iterator last,
				    Extractor extractor)
    {
      std::string rule = lhs;
      
      if (last != first) {
	rule += '(';
	for (/**/; first != last - 1; ++ first)
	  rule += extractor(*first) + '|';
	rule += extractor(*first) + ')';
      }
      return rule;
    }
    
    template <typename Iterator, typename Extractor>
    std::string extract_rule(const std::string& lhs,
			     Iterator first, Iterator last,
			     Extractor extractor)
    {
      std::string rule = lhs;
      
      if (last != first) {
	rule += '(';
	for (/**/; first != last - 1; ++ first)
	  rule += extractor(first->non_terminal()) + '|';
	rule += extractor(first->non_terminal()) + ')';
      }
      return rule;
    }

    
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
	  const std::string rule_string  = extract_phrase_rule(edge.rule->lhs, phrase.begin(), phrase.end(), __extractor_none());
	  
	  std::string rule_cluster = rule_string;
	  std::string rule_prefix  = rule_string;
	  std::string rule_suffix  = rule_string;
	  std::string rule_digits  = rule_string;
	  
	  if (cluster)
	    rule_cluster = extract_phrase_rule(edge.rule->lhs, phrase.begin(), phrase.end(), __extractor<cluster_type>(*cluster));
	  if (stemmer_prefix)
	    rule_prefix = extract_phrase_rule(edge.rule->lhs, phrase.begin(), phrase.end(), __extractor<stemmer_type>(*stemmer_prefix));
	  if (stemmer_suffix)
	    rule_suffix = extract_phrase_rule(edge.rule->lhs, phrase.begin(), phrase.end(), __extractor<stemmer_type>(*stemmer_suffix));
	  if (stemmer_digits)
	    rule_digits = extract_phrase_rule(edge.rule->lhs, phrase.begin(), phrase.end(), __extractor<stemmer_type>(*stemmer_digits));
	  
	  const id_type id_string  = rule_id(rule_string);
	  const id_type id_cluster = (cluster         ? rule_id(rule_cluster) : id_string);
	  const id_type id_prefix  = (stemmer_prefix  ? rule_id(rule_prefix)  : id_string);
	  const id_type id_suffix  = (stemmer_suffix  ? rule_id(rule_suffix)  : id_string);
	  const id_type id_digits  = (stemmer_digits  ? rule_id(rule_digits)  : id_string);
	  
	  apply_feature(features, rule_string);
	  if (id_cluster != id_string)
	    apply_feature(features, rule_cluster);
	  if (id_prefix != id_string)
	    apply_feature(features, rule_prefix);
	  if (id_suffix != id_string)
	    apply_feature(features, rule_suffix);
	  if (id_digits != id_string)
	    apply_feature(features, rule_digits);

	  id_type* context = reinterpret_cast<id_type*>(state);
	  context[0] = id_string;
	  context[1] = id_cluster;
	  context[2] = id_prefix;
	  context[3] = id_suffix;
	  context[4] = id_digits;
	} else {
	  int pos_non_terminal = 0;
	  phrase_type::const_iterator piter_end = phrase.end();
	  for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter)
	    if (piter->is_non_terminal()) {
	      int antecedent_index = piter->non_terminal_index() - 1;
	      if (antecedent_index < 0)
		antecedent_index = pos_non_terminal;
	      
	      const id_type* antecedent_context = reinterpret_cast<const id_type*>(states[antecedent_index]);
	      
	      apply_feature(features, edge.rule->lhs, string_map[antecedent_context[0]]);
	      if (antecedent_context[1] != antecedent_context[0])
		apply_feature(features, edge.rule->lhs, string_map[antecedent_context[1]]);
	      if (antecedent_context[2] != antecedent_context[0])
		apply_feature(features, edge.rule->lhs, string_map[antecedent_context[2]]);
	      if (antecedent_context[3] != antecedent_context[0])
		apply_feature(features, edge.rule->lhs, string_map[antecedent_context[3]]);
	      if (antecedent_context[4] != antecedent_context[0])
		apply_feature(features, edge.rule->lhs, string_map[antecedent_context[4]]);

	      ++ pos_non_terminal;
	    }
	  
	  const std::string rule_string  = extract_rule(edge.rule->lhs, phrase.begin(), phrase.end(), __extractor_none());
	  
	  std::string rule_cluster = rule_string;
	  std::string rule_prefix  = rule_string;
	  std::string rule_suffix  = rule_string;
	  std::string rule_digits  = rule_string;
	  
	  if (cluster)
	    rule_cluster = extract_rule(edge.rule->lhs, phrase.begin(), phrase.end(), __extractor<cluster_type>(*cluster));
	  if (stemmer_prefix)
	    rule_prefix = extract_rule(edge.rule->lhs, phrase.begin(), phrase.end(), __extractor<stemmer_type>(*stemmer_prefix));
	  if (stemmer_suffix)
	    rule_suffix = extract_rule(edge.rule->lhs, phrase.begin(), phrase.end(), __extractor<stemmer_type>(*stemmer_suffix));
	  if (stemmer_digits)
	    rule_digits = extract_rule(edge.rule->lhs, phrase.begin(), phrase.end(), __extractor<stemmer_type>(*stemmer_digits));
	  
	  const id_type id_string  = rule_id(rule_string);
	  const id_type id_cluster = (cluster         ? rule_id(rule_cluster) : id_string);
	  const id_type id_prefix  = (stemmer_prefix  ? rule_id(rule_prefix)  : id_string);
	  const id_type id_suffix  = (stemmer_suffix  ? rule_id(rule_suffix)  : id_string);
	  const id_type id_digits  = (stemmer_digits  ? rule_id(rule_digits)  : id_string);
	  
	  apply_feature(features, rule_string);
	  if (id_cluster != id_string)
	    apply_feature(features, rule_cluster);
	  if (id_prefix != id_string)
	    apply_feature(features, rule_prefix);
	  if (id_suffix != id_string)
	    apply_feature(features, rule_suffix);
	  if (id_digits != id_string)
	    apply_feature(features, rule_digits);

	  id_type* context = reinterpret_cast<id_type*>(state);
	  context[0] = id_string;
	  context[1] = id_cluster;
	  context[2] = id_prefix;
	  context[3] = id_suffix;
	  context[4] = id_digits;
	}
      }

      id_type rule_id(const std::string& rule) const
      {
	string_map_type::iterator iter = const_cast<string_map_type&>(string_map).insert(rule).first;
	
	return iter - string_map.begin();
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
      
      int stemmer_prefix_size = 0;
      int stemmer_suffix_size = 0;
      bool stemmer_digits = false;
      
      boost::filesystem::path cluster_path;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	  const std::string& yield = piter->second;
	  
	  if (strcasecmp(yield.c_str(), "source") == 0)
	    source = true;
	  else if (strcasecmp(yield.c_str(), "target") == 0)
	    target = true;
	  else
	    throw std::runtime_error("unknown parameter: " + parameter);
	} else if (strcasecmp(piter->first.c_str(), "cluster") == 0)
	  cluster_path = piter->second;
	else if (strcasecmp(piter->first.c_str(), "prefix") == 0)
	  stemmer_prefix_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "suffix") == 0)
	  stemmer_suffix_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "digits") == 0)
	  stemmer_digits = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for parent: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (source && target)
	throw std::runtime_error("both source and target?");
      if (! source && ! target)
	throw std::runtime_error("what side are you going to use?");
      
      if (stemmer_prefix_size < 0)
	throw std::runtime_error("negative prefix size?");
      if (stemmer_suffix_size < 0)
	throw std::runtime_error("negative suffix size?");

      std::auto_ptr<impl_type> parent_impl(source
					   ? dynamic_cast<impl_type*>(new __ParentImpl<__parent_extract_source>())
					   : dynamic_cast<impl_type*>(new __ParentImpl<__parent_extract_target>()));
      
      
      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file: " + cluster_path.file_string());
	
	parent_impl->cluster = &cicada::Cluster::create(cluster_path);
      }
      
      if (stemmer_prefix_size > 0)
	parent_impl->stemmer_prefix = &cicada::Stemmer::create("prefix:size=" + boost::lexical_cast<std::string>(stemmer_prefix_size));
      
      if (stemmer_suffix_size > 0)
	parent_impl->stemmer_suffix = &cicada::Stemmer::create("suffix:size=" + boost::lexical_cast<std::string>(stemmer_suffix_size));

      if (stemmer_digits)
	parent_impl->stemmer_digits = &cicada::Stemmer::create("digits");
      
      // parent conext (surface, cluster, prefix, suffix, digits)
      base_type::__state_size = sizeof(impl_type::id_type) * 5;
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
	pimpl = new parent_source_type(*dynamic_cast<const parent_source_type*>(x.pimpl));
      else
	pimpl = new parent_target_type(*dynamic_cast<const parent_target_type*>(x.pimpl));
    }
    
    Parent& Parent::operator=(const Parent& x)
    {
      typedef __ParentImpl<__parent_extract_source> parent_source_type;
      typedef __ParentImpl<__parent_extract_target> parent_target_type;

      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      std::auto_ptr<impl_type> tmp(pimpl);
      
      if (dynamic_cast<const parent_source_type*>(x.pimpl))
	pimpl = new parent_source_type(*dynamic_cast<const parent_source_type*>(x.pimpl));
      else
	pimpl = new parent_target_type(*dynamic_cast<const parent_target_type*>(x.pimpl));
      
      return *this;
    }


    template <typename FeaturePrefix, typename Feature>
    inline
    bool equal_prefix(const FeaturePrefix& prefix, const Feature& x)
    {
      return x.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), x.begin());
    }
    
    void Parent::apply(state_ptr_type& state,
		       const state_ptr_set_type& states,
		       const edge_type& edge,
		       feature_set_type& features,
		       feature_set_type& estimates,
		       const bool final) const
    {
      const std::string& __feature_prefix = base_type::feature_name();
      for (feature_set_type::iterator fiter = features.begin(); fiter != features.end(); /**/)
	if (equal_prefix(__feature_prefix, fiter->first))
	  features.erase(fiter ++);
	else
	  ++ fiter;
      
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();

      pimpl->parent_score(state, states, edge, features);

      if (final)
	pimpl->parent_final_score(state, edge, features);
    }

    void Parent::apply_coarse(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const
    {
    }

    void Parent::initialize()
    {
      pimpl->clear();
    }
  };
};
