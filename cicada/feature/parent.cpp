//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/parent.hpp"
#include "cicada/parameter.hpp"

#include "utils/indexed_set.hpp"
#include "utils/lexical_cast.hpp"

#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace feature
  {

    struct __extractor_none
    {
      template <typename Word>
      const std::string& operator()(const Word& word) const
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
	return extract(word);
      }
    };

    class ParentImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      typedef cicada::Cluster  cluster_type;
      typedef cicada::Stemmer  stemmer_type;

      typedef cicada::ClusterStemmer normalizer_type;
      typedef std::vector<normalizer_type, std::allocator<normalizer_type> > normalizer_set_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;
      
      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
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

      typedef utils::simple_vector<std::string, std::allocator<std::string> > normalized_set_type;
      typedef utils::chunk_vector<normalized_set_type, 4096 /sizeof(normalized_set_type), std::allocator<normalized_set_type> > normalized_map_type;
      
      ParentImpl()
	: exclude_terminal(false),
	  sentence(0),
	  forced_feature(false),
	  alignment_mode(false),
	  source_root_mode(false),
	  attr_target_position("target-position"),
	  attr_source_root("source-root") {}
      
      void clear()
      {
	string_map.clear();
	normalized.clear();
      }

      template <typename Iterator, typename Extractor>
      std::string extract_phrase_rule(const std::string& lhs,
				      Iterator first, Iterator last,
				      Extractor extractor) const
      {
	std::string rule = lhs;
      
	if (exclude_terminal) {
	  if (last != first) {
	    rule += '(';
	    for (/**/; first != last - 1; ++ first)
	      if (first->is_non_terminal())
		rule += extractor(*first) + '|';
	    if (first->is_non_terminal())
	      rule += extractor(*first);
	    rule += ')';
	  }
	} else {
	  if (last != first) {
	    rule += '(';
	    for (/**/; first != last - 1; ++ first)
	      rule += extractor(*first) + '|';
	    rule += extractor(*first) + ')';
	  }
	}

	return rule;
      }
    
      template <typename Iterator, typename Extractor>
      std::string extract_rule(const std::string& lhs,
			       Iterator first, Iterator last,
			       Extractor extractor) const
      {
	std::string rule = lhs;

	if (exclude_terminal) {
	  if (last != first) {
	    rule += '(';
	    for (/**/; first != last - 1; ++ first)
	      if (first->is_non_terminal())
		rule += extractor(first->non_terminal()) + '|';
	    if (first->is_non_terminal())
	      rule += extractor(first->non_terminal());
	    rule += ')';
	  }
	} else {
	  if (last != first) {
	    rule += '(';
	    for (/**/; first != last - 1; ++ first)
	      rule += extractor(first->non_terminal()) + '|';
	    rule += extractor(first->non_terminal()) + ')';
	  }
	}
	return rule;
      }

      normalizer_set_type normalizers;

      bool exclude_terminal;
      
      string_map_type     string_map;
      normalized_map_type normalized;

      feature_type feature_name_prefix;

      const sentence_type* sentence;
      
      bool forced_feature;
      bool alignment_mode;
      bool source_root_mode;
      
      attribute_type attr_target_position;
      attribute_type attr_source_root;

      sentence_type __buffer;

      struct __attribute_integer : public boost::static_visitor<cicada::AttributeVector::int_type>
      {
	typedef cicada::AttributeVector attribute_set_type;
	
	attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
	attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -2; }
	attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -2; }
      };
      
      struct __attribute_string : public boost::static_visitor<cicada::AttributeVector::string_type>
      {
	typedef cicada::AttributeVector attribute_set_type;
	
	attribute_set_type::string_type operator()(const attribute_set_type::int_type& x) const { return ""; }
	attribute_set_type::string_type operator()(const attribute_set_type::float_type& x) const { return ""; }
	attribute_set_type::string_type operator()(const attribute_set_type::string_type& x) const { return x; }
      };

      symbol_type root_label(const std::string& rule) const
      {
	std::string::size_type pos = rule.find('(');
	if (pos != std::string::npos) {
	  const std::string label = rule.substr(0, pos);
	  if (label.empty() || label[0] != '[' || label[label.size() -1] != ']')
	    throw std::runtime_error("invlaid label:" + label);
	  return label;
	} else {
	  if (rule.empty() || rule[0] != '[' || rule[rule.size() -1] != ']')
	    throw std::runtime_error("invlaid label: " + rule);
	  
	  return rule;
	}
      }
      
      symbol_type root_label(const edge_type& edge) const
      {
	if (source_root_mode) {
	  std::string label;
	  
	  attribute_set_type::const_iterator riter = edge.attributes.find(attr_source_root);
	  if (riter != edge.attributes.end())
	    label = boost::apply_visitor(__attribute_string(), riter->second);
	  
	  if (label.empty())
	    return edge.rule->lhs;
	  else
	    return label;
	} else
	  return edge.rule->lhs;
      }
      
      void parent_score(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features) const
      {
	// this feature function is complicated in that we know nothing about the source-side...
	
	const rule_type::symbol_set_type& phrase = edge.rule->rhs;
	
	if (states.empty()) {
	  const symbol_type cat = root_label(edge);
	  
	  if (alignment_mode) {
	    symbol_type symbol = vocab_type::EPSILON;
	    
	    attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_position);
	    if (titer == edge.attributes.end())
	      throw std::runtime_error("we do not support non alignment forest");
	    
	    const int target_pos = boost::apply_visitor(__attribute_integer(), titer->second);
	    if (sentence && target_pos >= 0)
	      symbol = sentence->operator[](target_pos);
	    
	    const std::string rule_string  = extract_phrase_rule(cat, &symbol, (&symbol) + 1, __extractor_none());
	    const id_type     id_string    = rule_id(rule_string);
	    
	    apply_feature(features, rule_string);
	    
	    for (size_t i = 0; i != normalizers.size(); ++ i) {
	      std::string& rule_norm = const_cast<std::string&>(normalized[id_string][i]);
	      
	      if (rule_norm.empty())
		rule_norm = extract_phrase_rule(cat, &symbol, (&symbol) + 1, __extractor<normalizer_type>(normalizers[i]));
	      
	      if (rule_string != rule_norm)
		apply_feature(features, rule_norm);
	    }
	    
	    *reinterpret_cast<id_type*>(state) = id_string;
	  } else {
	    const std::string rule_string  = extract_phrase_rule(cat, phrase.begin(), phrase.end(), __extractor_none());
	    const id_type     id_string    = rule_id(rule_string);
	    
	    apply_feature(features, rule_string);
	    
	    for (size_t i = 0; i != normalizers.size(); ++ i) {
	      std::string& rule_norm = const_cast<std::string&>(normalized[id_string][i]);
	      
	      if (rule_norm.empty())
		rule_norm = extract_phrase_rule(cat, phrase.begin(), phrase.end(), __extractor<normalizer_type>(normalizers[i]));
	      
	      if (rule_string != rule_norm)
		apply_feature(features, rule_norm);
	    }
	    
	    *reinterpret_cast<id_type*>(state) = id_string;
	  }
	} else {
	  const symbol_type cat = root_label(edge);
	  
	  sentence_type& buffer = const_cast<sentence_type&>(__buffer);
	  buffer.clear();
	  
	  int pos_non_terminal = 0;
	  phrase_type::const_iterator piter_end = phrase.end();
	  for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter)
	    if (piter->is_non_terminal()) {
	      int antecedent_index = piter->non_terminal_index() - 1;
	      if (antecedent_index < 0)
		antecedent_index = pos_non_terminal;
	      ++ pos_non_terminal;
	      
	      const id_type* antecedent_context = reinterpret_cast<const id_type*>(states[antecedent_index]);
	      
	      buffer.push_back(root_label(string_map[*antecedent_context]));
	      
	      apply_feature(features, cat, string_map[*antecedent_context]);
	      
	      for (size_t i = 0; i != normalizers.size(); ++ i)
		if (string_map[*antecedent_context] != normalized[*antecedent_context][i])
		  apply_feature(features, cat, normalized[*antecedent_context][i]);
	    } else
	      buffer.push_back(*piter);
	  
	  const std::string rule_string  = extract_rule(cat, buffer.begin(), buffer.end(), __extractor_none());
	  const id_type id_string        = rule_id(rule_string);
	  
	  apply_feature(features, rule_string);
	  
	  for (size_t i = 0; i != normalizers.size(); ++ i) {
	    std::string& rule_norm = const_cast<std::string&>(normalized[id_string][i]);
	    
	    if (rule_norm.empty())
	      rule_norm = extract_rule(cat, buffer.begin(), buffer.end(), __extractor<normalizer_type>(normalizers[i]));
	    
	    if (rule_string != rule_norm)
	      apply_feature(features, rule_norm);
	  }
	  
	  *reinterpret_cast<id_type*>(state) = id_string;
	}
      }

      id_type rule_id(const std::string& rule) const
      {
	string_map_type::iterator iter = const_cast<string_map_type&>(string_map).insert(rule).first;

	const id_type id = iter - string_map.begin();

	if (! normalizers.empty()) {
	  if (id >= normalized.size())
	    const_cast<normalized_map_type&>(normalized).resize(id + 1);

	  if (normalized[id].empty())
	    const_cast<normalized_set_type&>(normalized[id]) = normalized_set_type(normalizers.size());
	}
	
	return id;
      }
      
      
      void parent_final_score(const state_ptr_type& state,
			      const edge_type& edge,
			      feature_set_type& features) const
      {
	const id_type* id = reinterpret_cast<const id_type*>(state);
	
	apply_feature(features, vocab_type::BOS, string_map[*id]);
      }
      
      void apply_feature(feature_set_type& features,
			 const std::string& rule) const
      {
	const std::string name = static_cast<const std::string&>(feature_name_prefix) + ":" + rule;
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;
      }
      
      void apply_feature(feature_set_type& features,
			 const std::string& parent,
			 const std::string& rule) const
      {
	const std::string name = static_cast<const std::string&>(feature_name_prefix) + ":" + parent + '(' + rule + ')';
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;
      }
    };
    
    Parent::Parent(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "parent")
	throw std::runtime_error("is this really parent feature function? " + parameter);

      impl_type::normalizer_set_type normalizers;
      std::string name;
      bool alignment_mode = false;
      bool source_root_mode = false;
      bool exclude_terminal = false;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "cluster") == 0) {
	  if (! boost::filesystem::exists(piter->second))
	    throw std::runtime_error("no cluster file: " + piter->second);
	  
	  normalizers.push_back(impl_type::normalizer_type(&cicada::Cluster::create(piter->second)));
	} else if (strcasecmp(piter->first.c_str(), "stemmer") == 0)
	  normalizers.push_back(impl_type::normalizer_type(&cicada::Stemmer::create(piter->second)));
	else if (strcasecmp(piter->first.c_str(), "exclude-terminal") == 0)
	  exclude_terminal = utils::lexical_cast<bool>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "name") == 0)
	  name = piter->second;
	else if (strcasecmp(piter->first.c_str(), "alignment") == 0)
	  alignment_mode = utils::lexical_cast<bool>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "source-root") == 0)
	  source_root_mode = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for parent: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> parent_impl(new impl_type());
      
      parent_impl->normalizers.swap(normalizers);
      
      parent_impl->exclude_terminal = exclude_terminal;
      parent_impl->alignment_mode = alignment_mode;
      parent_impl->source_root_mode = source_root_mode;
      parent_impl->feature_name_prefix = (name.empty() ? std::string("parent") : name);
      
      // parent conext 
      base_type::__state_size = sizeof(impl_type::id_type);
      base_type::__feature_name = (name.empty() ? std::string("parent") : name);
      base_type::__sparse_feature = true;
      
      pimpl = parent_impl.release();
    }
    
    Parent::~Parent() { std::auto_ptr<impl_type> tmp(pimpl); }

    
    Parent::Parent(const Parent& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}
    
    Parent& Parent::operator=(const Parent& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      return *this;
    }
    
    void Parent::apply(state_ptr_type& state,
		       const state_ptr_set_type& states,
		       const edge_type& edge,
		       feature_set_type& features,
		       feature_set_type& estimates,
		       const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
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
    void Parent::apply_predict(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       feature_set_type& estimates,
			       const bool final) const
    {}
    void Parent::apply_scan(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    const int dot,
			    feature_set_type& features,
			    feature_set_type& estimates,
			    const bool final) const
    {}
    void Parent::apply_complete(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }


    void Parent::initialize()
    {
      pimpl->clear();
    }

    void Parent::assign(const size_type& id,
			const hypergraph_type& hypergraph,
			const lattice_type& lattice,
			const span_set_type& spans,
			const sentence_set_type& targets,
			const ngram_count_set_type& ngram_counts)
    {
      pimpl->sentence = 0;
      if (! targets.empty())
	pimpl->sentence = &targets.front();
    }
    
  };
};
