//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/antecedent.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"

#include "utils/compact_trie_dense.hpp"
#include "utils/lexical_cast.hpp"

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace feature
  {
    class AntecedentImpl
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
      
      typedef utils::compact_trie_dense<symbol_type, std::string, boost::hash<symbol_type>, std::equal_to<symbol_type>,
					std::allocator<std::pair<const symbol_type, std::string> > > tree_map_type;
      
      typedef tree_map_type::id_type id_type;
      
      AntecedentImpl()
	: tree_map(symbol_type()),
	  sentence(0),
	  forced_feature(false),
	  alignment_mode(false),
	  attr_target_position("target-position") {}
      
      void clear()
      {
	tree_map.clear();
      }

      normalizer_set_type normalizers;
      
      tree_map_type  tree_map;

      feature_type feature_name_prefix;

      const sentence_type* sentence;

      bool forced_feature;
      bool alignment_mode;

      attribute_type attr_target_position;
      
      struct __attribute_integer : public boost::static_visitor<cicada::AttributeVector::int_type>
      {
	typedef cicada::AttributeVector attribute_set_type;
	
	attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
	attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -2; }
	attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -2; }
      };

      void antecedent_score(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features) const
      {
	// this feature function is complicated in that we know nothing about the source-side...
	
	const rule_type::symbol_set_type& phrase = edge.rule->rhs;
	
	if (states.empty()) {
	  // we do not add feature here, since we know nothing abount surrounding context...
	  symbol_type prefix = vocab_type::EPSILON;
	  symbol_type suffix = vocab_type::EPSILON;
	  int span_size = 0;
	  
	  if (alignment_mode) {
	    attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_position);
	    if (titer == edge.attributes.end())
	      throw std::runtime_error("we do not support non alignment forest");
	    
	    const int target_pos = boost::apply_visitor(__attribute_integer(), titer->second);
	    
	    if (sentence && target_pos >= 0) {
	      const symbol_type& target = sentence->operator[](target_pos);
	      
	      prefix = target;
	      suffix = target;
	      span_size = 1;
	    }
	    
	  } else {
	    phrase_type::const_iterator piter_end = phrase.end();
	    for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter)
	      if (*piter != vocab_type::EPSILON) {
		if (prefix == vocab_type::EPSILON)
		  prefix = *piter;
		suffix = *piter;
		++ span_size;
	      }
	  }
	  
	  id_type*     context_tree   = reinterpret_cast<id_type*>(state);
	  symbol_type* context_symbol = reinterpret_cast<symbol_type*>(context_tree + 1);
	  int*         context_size   = reinterpret_cast<int*>(context_symbol + 2);
	  
	  *context_tree = tree_map.root();
	  context_symbol[0] = prefix;
	  context_symbol[1] = suffix;
	  *context_size = span_size;
	} else {
	  symbol_type prefix = vocab_type::EMPTY;
	  symbol_type suffix = vocab_type::EMPTY;
	  
	  std::string antecedent_string;
	  id_type     node = tree_map.root();
	  int span_size = 0;
	  
	  int pos_non_terminal = 0;
	  phrase_type::const_iterator piter_end = phrase.end();
	  for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter)
	    if (piter->is_non_terminal()) {
	      int antecedent_index = piter->non_terminal_index() - 1;
	      if (antecedent_index < 0)
		antecedent_index = pos_non_terminal;
	      
	      const id_type*     antecedent_tree   = reinterpret_cast<const id_type*>(states[antecedent_index]);
	      const symbol_type* antecedent_symbol = reinterpret_cast<const symbol_type*>(antecedent_tree + 1);
	      const int*         antecedent_size   = reinterpret_cast<const int*>(antecedent_symbol + 2);
	      
	      node = tree_id(*piter, node);
	      
	      antecedent_string += compose_tree(*piter, *antecedent_tree);
	      span_size += *antecedent_size;
	      
	      if (prefix == vocab_type::EMPTY)
		prefix = antecedent_symbol[0];
	      suffix = antecedent_symbol[1];
	      
	      ++ pos_non_terminal;
	    } else if (*piter != vocab_type::EPSILON) {
	      ++ span_size;
	      
	      if (prefix == vocab_type::EMPTY)
		prefix = *piter;
	      suffix = *piter;
	    }
	  
	  // apply feature...
	  apply_feature(features, edge.rule->lhs, antecedent_string, prefix, suffix, span_size);
	  
	  // next context...
	  id_type*     context_tree   = reinterpret_cast<id_type*>(state);
	  symbol_type* context_symbol = reinterpret_cast<symbol_type*>(context_tree + 1);
	  int*         context_size   = reinterpret_cast<int*>(context_symbol + 2);
	  
	  *context_tree = node;
	  context_symbol[0] = prefix;
	  context_symbol[1] = suffix;
	  *context_size = span_size;
	}
      }
      
      void antecedent_final_score(const state_ptr_type& state,
				  feature_set_type& features) const
      {
	// nothing to apply!
      }

      const std::string compose_tree(const std::string& node, const id_type& id) const
      {
	if (tree_map.is_root(id))
	  return '(' + node + ')';
	else
	  return '(' + node + '(' + tree_map[id] + "))";
      }

      
      id_type tree_id(const symbol_type& node, const id_type parent) const
      {
	tree_map_type& __tree_map = const_cast<tree_map_type&>(tree_map);
	
	const id_type id = __tree_map.insert(parent, node);
	
	if (__tree_map[id].empty()) {
	  if (! __tree_map.is_root(parent))
	    __tree_map[id] =  __tree_map[parent] + static_cast<const std::string&>(node);
	  else
	    __tree_map[id] = node;
	}
	
	return id;
      }

      void apply_feature(feature_set_type& features,
			 const std::string& node,
			 const std::string& antecedent,
			 const symbol_type& prefix, const symbol_type& suffix,
			 const int span_size) const
      {
	const std::string name = feature_name(node, antecedent, prefix, suffix, span_size);
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;
	
	for (size_t i = 0; i != normalizers.size(); ++ i) {
	  const symbol_type prefix_norm = normalizers[i](prefix);
	  const symbol_type suffix_norm = normalizers[i](suffix);
	  
	  if (prefix_norm != prefix || suffix_norm != suffix) {
	    const std::string name = feature_name(node, antecedent, prefix_norm, suffix_norm, span_size);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}
      }
      
      const std::string feature_name(const std::string& node,
				     const std::string& antecedent,
				     const std::string& prefix,
				     const std::string& suffix,
				     const int span_size) const
      {
	return (static_cast<const std::string&>(feature_name_prefix) + ":"
		+ node + antecedent + '|' + prefix + '|' + suffix + '|' + boost::lexical_cast<std::string>(span_size));
      }
    };
    
    
    Antecedent::Antecedent(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "antecedent")
	throw std::runtime_error("is this really antecedent feature function? " + parameter);

      impl_type::normalizer_set_type normalizers;
      std::string name;
      bool alignment_mode = false;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "cluster") == 0) {
	  if (! boost::filesystem::exists(piter->second))
	    throw std::runtime_error("no cluster file: " + piter->second);
	  
	  normalizers.push_back(impl_type::normalizer_type(&cicada::Cluster::create(piter->second)));
	} else if (strcasecmp(piter->first.c_str(), "stemmer") == 0)
	  normalizers.push_back(impl_type::normalizer_type(&cicada::Stemmer::create(piter->second)));
	else if (strcasecmp(piter->first.c_str(), "name") == 0)
	  name = piter->second;
	else if (strcasecmp(piter->first.c_str(), "alignment") == 0)
	  alignment_mode = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for antecedent: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> antecedent_impl(new impl_type());
      
      antecedent_impl->normalizers.swap(normalizers);
      antecedent_impl->alignment_mode = alignment_mode;
      antecedent_impl->feature_name_prefix = (name.empty() ? std::string("antecedent") : name);
      
      // antecedent conext + terminal-boundary + span-size
      base_type::__state_size = sizeof(impl_type::id_type) + sizeof(symbol_type) * 2 + sizeof(int);
      base_type::__feature_name = (name.empty() ? std::string("antecedent") : name);
      base_type::__sparse_feature = true;
      
      pimpl = antecedent_impl.release();
    }
    
    Antecedent::~Antecedent() { std::auto_ptr<impl_type> tmp(pimpl); }

    
    Antecedent::Antecedent(const Antecedent& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    { }
    
    
    Antecedent& Antecedent::operator=(const Antecedent& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void Antecedent::apply(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge,
			   feature_set_type& features,
			   feature_set_type& estimates,
			   const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();

      pimpl->antecedent_score(state, states, edge, features);
      
      if (final)
	pimpl->antecedent_final_score(state, features);
    }

    void Antecedent::apply_coarse(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {
      
    }
    
    void Antecedent::apply_predict(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {}
    
    void Antecedent::apply_scan(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				const int dot,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {}
    void Antecedent::apply_complete(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    feature_set_type& estimates,
				    const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    void Antecedent::initialize()
    {
      pimpl->clear();
    }

    void Antecedent::assign(const size_type& id,
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
