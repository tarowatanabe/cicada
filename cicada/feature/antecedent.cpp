//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/antecedent.hpp"
#include "cicada/feature/feature_builder.hpp"

#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"

#include "utils/trie_compact.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/piece.hpp"
#include "utils/bithack.hpp"

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
      
      typedef utils::trie_compact<symbol_type, std::string,
				  utils::unassigned<symbol_type>, 
				  boost::hash<symbol_type>, std::equal_to<symbol_type>,
				  std::allocator<std::pair<const symbol_type, std::string> > > tree_map_type;
      
      typedef tree_map_type::id_type id_type;

      typedef FeatureBuilder feature_builder_type;
      
      AntecedentImpl()
	: tree_map(),
	  sentence(0),
	  forced_feature(false),
	  alignment_mode(false),
	  source_root_mode(false),
	  attr_target_position("target-position"),
	  attr_source_root("source-root") {}
      
      void clear()
      {
	tree_map.clear();
      }

      normalizer_set_type normalizers;
      
      tree_map_type  tree_map;

      feature_type feature_name_prefix;
      
      feature_builder_type feature_builder;
      feature_builder_type tree_builder;

      const sentence_type* sentence;

      bool forced_feature;
      bool alignment_mode;
      bool source_root_mode;

      attribute_type attr_target_position;
      attribute_type attr_source_root;
      
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
	  symbol_type* context_root   = reinterpret_cast<symbol_type*>(context_size + 1);
	  
	  *context_tree = tree_map.root();
	  context_symbol[0] = prefix;
	  context_symbol[1] = suffix;
	  *context_size = span_size;
	  *context_root = root_label(edge);
	} else {
	  symbol_type prefix = vocab_type::EMPTY;
	  symbol_type suffix = vocab_type::EMPTY;
	  
	  feature_builder_type& builder = const_cast<feature_builder_type&>(tree_builder);
	  builder.clear();
	  
	  //std::string antecedent_string;
	  id_type node = tree_map.root();
	  int     span_size = 0;
	  
	  int pos_non_terminal = 0;
	  phrase_type::const_iterator piter_end = phrase.end();
	  for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter)
	    if (piter->is_non_terminal()) {
	      const int __non_terminal_index = piter->non_terminal_index();
	      const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, pos_non_terminal, __non_terminal_index - 1);
	      
	      const id_type*     antecedent_tree   = reinterpret_cast<const id_type*>(states[antecedent_index]);
	      const symbol_type* antecedent_symbol = reinterpret_cast<const symbol_type*>(antecedent_tree + 1);
	      const int*         antecedent_size   = reinterpret_cast<const int*>(antecedent_symbol + 2);
	      const symbol_type* antecedent_root   = reinterpret_cast<const symbol_type*>(antecedent_size + 1);
	      
	      node = tree_id(*antecedent_root, node);
	      //antecedent_string += compose_tree(*antecedent_root, *antecedent_tree);
	      
	      if (tree_map.is_root(*antecedent_tree))
		builder << "(" << *antecedent_root << ")";
	      else
		builder << "(" << *antecedent_root << "(" << tree_map[*antecedent_tree] << "))";

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
	  const symbol_type cat = root_label(edge);

	  //apply_feature(features, cat, antecedent_string, prefix, suffix, span_size);
	  apply_feature(features, cat, builder, prefix, suffix, span_size);
	  
	  // next context...
	  id_type*     context_tree   = reinterpret_cast<id_type*>(state);
	  symbol_type* context_symbol = reinterpret_cast<symbol_type*>(context_tree + 1);
	  int*         context_size   = reinterpret_cast<int*>(context_symbol + 2);
	  symbol_type* context_root   = reinterpret_cast<symbol_type*>(context_size + 1);
	  
	  *context_tree = node;
	  context_symbol[0] = prefix;
	  context_symbol[1] = suffix;
	  *context_size = span_size;
	  *context_root = cat;
	}
      }
      
      void antecedent_final_score(const state_ptr_type& state,
				  feature_set_type& features) const
      {
	// nothing to apply!
      }

#if 0
      const std::string compose_tree(const std::string& node, const id_type& id) const
      {
	if (tree_map.is_root(id))
          return '(' + node + ')';
        else
          return '(' + node + '(' + tree_map[id] + "))";
      }
#endif
      
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
			 const feature_builder_type& antecedent,
			 const symbol_type& prefix, const symbol_type& suffix,
			 const int span_size) const
      {
	const int span_size_power2 = utils::bithack::branch(utils::bithack::is_power2(span_size),
							    span_size,
							    static_cast<int>(utils::bithack::next_largest_power2(span_size)));
	
	feature_builder_type& builder = const_cast<feature_builder_type&>(feature_builder);
	builder.clear();
	
	builder << feature_name_prefix
		<< ":" << node << antecedent
		<< "|" << prefix
		<< "|" << suffix
		<< "|" << span_size_power2;

	if (forced_feature || builder.exists())
	  features[builder] += 1.0;
	
	for (size_t i = 0; i != normalizers.size(); ++ i) {
	  builder.clear();
	  builder << feature_name_prefix
		  << ":" << node << antecedent
		  << "|" << normalizers[i](prefix)
		  << "|" << normalizers[i](suffix)
		  << "|" << span_size_power2;
	  
	  if (forced_feature || builder.exists())
	    features[builder] += 1.0;
	}
      }
    };
    
    
    Antecedent::Antecedent(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (utils::ipiece(param.name()) != "antecedent")
	throw std::runtime_error("is this really antecedent feature function? " + parameter);

      impl_type::normalizer_set_type normalizers;
      std::string name;
      bool alignment_mode = false;
      bool source_root_mode = false;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "cluster") {
	  if (! boost::filesystem::exists(piter->second))
	    throw std::runtime_error("no cluster file: " + piter->second);
	  
	  normalizers.push_back(impl_type::normalizer_type(&cicada::Cluster::create(piter->second)));
	} else if (utils::ipiece(piter->first) == "stemmer")
	  normalizers.push_back(impl_type::normalizer_type(&cicada::Stemmer::create(piter->second)));
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else if (utils::ipiece(piter->first) == "alignment")
	  alignment_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "source-root")
	  source_root_mode = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for antecedent: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::unique_ptr<impl_type> antecedent_impl(new impl_type());
      
      antecedent_impl->normalizers.swap(normalizers);
      antecedent_impl->alignment_mode = alignment_mode;
      antecedent_impl->source_root_mode = source_root_mode;
      antecedent_impl->feature_name_prefix = (name.empty() ? std::string("antecedent") : name);
      
      // antecedent conext + terminal-boundary + span-size
      base_type::__state_size = sizeof(impl_type::id_type) + sizeof(symbol_type) * 2 + sizeof(int) + sizeof(symbol_type);
      base_type::__feature_name = (name.empty() ? std::string("antecedent") : name);
      base_type::__sparse_feature = true;
      
      pimpl = antecedent_impl.release();
    }
    
    Antecedent::~Antecedent() { std::unique_ptr<impl_type> tmp(pimpl); }

    
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
			   const bool final) const
    {
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();

      feature_set_type feats;
      
      pimpl->antecedent_score(state, states, edge, feats);
      
      if (final)
	pimpl->antecedent_final_score(state, feats);

      features.update(feats, static_cast<const std::string&>(base_type::feature_name()));
    }

    void Antecedent::apply_coarse(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  const bool final) const
    {
      
    }
    
    void Antecedent::apply_predict(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   const bool final) const
    {}
    
    void Antecedent::apply_scan(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				const int dot,
				feature_set_type& features,
				const bool final) const
    {}
    void Antecedent::apply_complete(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    const bool final) const
    {
      apply(state, states, edge, features, final);
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
