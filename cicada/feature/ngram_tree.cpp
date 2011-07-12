//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/ngram_tree.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"

#include "utils/compact_trie_dense.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/piece.hpp"

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
      
      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;
      
      struct node_pair_type
      {
	typedef utils::simple_vector<std::string, std::allocator<std::string> > node_set_type;

	node_set_type nodes;
	
	node_pair_type() : nodes() {}
      };
      
      typedef utils::compact_trie_dense<symbol_type, node_pair_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
					std::allocator<std::pair<const symbol_type, node_pair_type> > > tree_map_type;

      typedef tree_map_type::id_type id_type;

      NGramTreeImpl()
	: tree_map(symbol_type()),
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
      
      phrase_span_set_type phrase_spans_impl;
      
      feature_type feature_name_prefix;
      
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
      
      void ngram_tree_score(state_ptr_type& state,
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

	  if (alignment_mode) {
	    attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_position);
	    if (titer == edge.attributes.end())
	      throw std::runtime_error("we do not support non alignment forest");
	    
	    const int target_pos = boost::apply_visitor(__attribute_integer(), titer->second);
	    
	    if (sentence && target_pos >= 0) {
	      const symbol_type& target = sentence->operator[](target_pos);
	      
	      compute_bound(&target, (&target) + 1, prefix, suffix);
	    }
	      
	  } else
	    compute_bound(phrase.begin(), phrase.end(), prefix, suffix);

	  const symbol_type cat = root_label(edge);
	  
	  id_type* context = reinterpret_cast<id_type*>(state);
	  context[0] = tree_id(cat, tree_id(prefix, tree_map.root()));
	  context[1] = tree_id(cat, tree_id(suffix, tree_map.root()));

	} else {
	  phrase_span_set_type& phrase_spans = const_cast<phrase_span_set_type&>(phrase_spans_impl);
	  
	  phrase_spans.clear();
	  phrase.terminals(std::back_inserter(phrase_spans));

	  if (phrase_spans.size() != states.size() + 1)
	    throw std::runtime_error("# of states does not match...");
	  
	  symbol_type prefix;
	  symbol_type suffix;
	  
	  compute_bound(phrase_spans.front().first, phrase_spans.front().second, prefix, suffix);
	  
	  id_type prefix_id = (prefix.empty() ? tree_map.root() : tree_id(prefix, tree_map.root()));
	  id_type suffix_id = (suffix.empty() ? tree_map.root() : tree_id(suffix, tree_map.root()));

	  const symbol_type cat = root_label(edge);
	  
	  phrase_span_set_type::const_iterator siter_begin = phrase_spans.begin();
	  phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	  for (phrase_span_set_type::const_iterator siter = siter_begin + 1; siter != siter_end; ++ siter) {
	    const phrase_span_type& span = *siter;
	    
	    // incase, we are working with non-synchronous parsing!
	    const int __non_terminal_index = (span.first - 1)->non_terminal_index();
	    const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, int(siter - (siter_begin + 1)), __non_terminal_index - 1);
	    
	    const id_type* antecedent_context = reinterpret_cast<const id_type*>(states[antecedent_index]);
	    //const symbol_type* antecedent_root = reinterpret_cast<const symbol_type*>(antecedent_context + 2);
	    
	    const id_type prefix_antecedent_id = antecedent_context[0];
	    const id_type suffix_antecedent_id = antecedent_context[1];
	    //const id_type prefix_antecedent_id = tree_id(*(span.first - 1), antecedent_context[0]);
	    //const id_type suffix_antecedent_id = tree_id(*(span.first - 1), antecedent_context[1]);
	    //const id_type prefix_antecedent_id = tree_id(*antecedent_root, antecedent_context[0]);
	    //const id_type suffix_antecedent_id = tree_id(*antecedent_root, antecedent_context[1]);
	    
	    symbol_type prefix_next;
	    symbol_type suffix_next;
	    compute_bound(span.first, span.second, prefix_next, suffix_next);
	    
	    const id_type prefix_next_id = (prefix_next.empty() ? tree_map.root() : tree_id(prefix_next, tree_map.root()));
	    const id_type suffix_next_id = (suffix_next.empty() ? tree_map.root() : tree_id(suffix_next, tree_map.root()));
	    
	    if (! tree_map.is_root(suffix_id))
	      apply_feature(features, cat, suffix_id, prefix_antecedent_id);
	    
	    if (tree_map.is_root(prefix_id))
	      prefix_id = prefix_antecedent_id;
	    
	    if (! tree_map.is_root(prefix_next_id)) {
	      apply_feature(features, cat, suffix_antecedent_id, prefix_next_id);
	      suffix_id = suffix_next_id;
	    } else
	      suffix_id = suffix_antecedent_id;
	  }
	  
	  // construct state...
	  id_type* context = reinterpret_cast<id_type*>(state);
	  
	  if (tree_map.is_root(prefix_id))
	    prefix_id = tree_id(vocab_type::EPSILON, tree_map.root());
	  if (tree_map.is_root(suffix_id))
	    suffix_id = tree_id(vocab_type::EPSILON, tree_map.root());

	  context[0] = tree_id(cat, prefix_id);
	  context[1] = tree_id(cat, suffix_id);
	}
      }

      void ngram_tree_final_score(const state_ptr_type& state,
				  const edge_type& edge,
				  feature_set_type& features) const
      {
	const id_type* antecedent_context = reinterpret_cast<const id_type*>(state);
	
	const id_type prefix_antecedent_id = antecedent_context[0];
	const id_type suffix_antecedent_id = antecedent_context[1];
	
	const symbol_type cat = root_label(edge);
	
	apply_feature(features, cat, tree_id(vocab_type::BOS, tree_map.root()), prefix_antecedent_id);
	apply_feature(features, cat, suffix_antecedent_id, tree_id(vocab_type::EOS, tree_map.root()));
      }

      
      id_type tree_id(const symbol_type& node, const id_type parent) const
      {
	tree_map_type& __tree_map = const_cast<tree_map_type&>(tree_map);
	
	const id_type id = __tree_map.insert(parent, node);
	
	if (__tree_map[id].nodes.empty())
	  __tree_map[id].nodes = node_pair_type::node_set_type(normalizers.size() + 1);
	
	if (__tree_map[id].nodes.front().empty()) {
	  if (__tree_map.is_root(parent)) {
	    __tree_map[id].nodes.front() = node;
	    
	    for (size_t i = 0; i != normalizers.size(); ++ i) 
	      __tree_map[id].nodes[i + 1] = normalizers[i](node);
	    
	  } else {
	    __tree_map[id].nodes.front() = compose_path(node, __tree_map[parent].nodes.front());
	    
	    for (size_t i = 0; i != normalizers.size(); ++ i) 
	      __tree_map[id].nodes[i + 1] = compose_path(node, __tree_map[parent].nodes[i + 1]);
	  }
	}
	
	return id;
      }

      void apply_feature(feature_set_type& features, const symbol_type& node, const id_type& prev, const id_type& next) const
      {
	const node_pair_type& prev_node = tree_map[prev];
	const node_pair_type& next_node = tree_map[next];
	
	const std::string name = feature_name(node, prev_node.nodes.front(), next_node.nodes.front());
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;

	for (size_t i = 0; i != normalizers.size(); ++ i) 
	  if (prev_node.nodes.front() != prev_node.nodes[i + 1] || next_node.nodes.front() != next_node.nodes[i + 1]) {
	    const std::string name = feature_name(node, prev_node.nodes[i + 1], next_node.nodes[i + 1]);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
      }

      template <typename Iterator>
      void compute_bound(Iterator first, Iterator last, symbol_type& prefix, symbol_type& suffix) const
      {
	for (Iterator iter = first; iter != last; ++ iter)
	  if (*iter != vocab_type::EPSILON) {
	    prefix = *iter;
	    break;
	  }
	
	for (Iterator iter = last; iter != first; -- iter)
	  if (*(iter - 1) != vocab_type::EPSILON) {
	    suffix = *(iter - 1);
	    break;
	  }
      }
      
      
      const std::string compose_path(const std::string& node, const std::string& antecedent) const
      {
	return node + '(' + antecedent + ')';
      }

      const std::string compose_tree(const std::string& node, const std::string& prev, const std::string& next) const
      {
	return node + '(' + prev + ")(" + next + ')';
      }
      
      const std::string feature_name(const std::string& node, const std::string& prev, const std::string& next) const
      {
	return static_cast<const std::string&>(feature_name_prefix) + ":" +  compose_tree(node, prev, next);
      }

      
    };

    
    NGramTree::NGramTree(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (utils::ipiece(param.name()) != "ngram-tree")
	throw std::runtime_error("is this really ngram tree feature function? " + parameter);
      
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
	  std::cerr << "WARNING: unsupported parameter for ngram-tree: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> ngram_tree_impl(new impl_type());

      ngram_tree_impl->normalizers.swap(normalizers);
      ngram_tree_impl->alignment_mode = alignment_mode;
      ngram_tree_impl->source_root_mode = source_root_mode;
      ngram_tree_impl->feature_name_prefix = (name.empty() ? std::string("ngram-tree") : name);
      
      // non-terminal + two neighbouring symbols + span-size
      base_type::__state_size = sizeof(impl_type::id_type) * 2;
      base_type::__feature_name = (name.empty() ? std::string("ngram-tree") : name);
      base_type::__sparse_feature = true;
      
      pimpl = ngram_tree_impl.release();
    }
    
    NGramTree::~NGramTree() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    
    NGramTree::NGramTree(const NGramTree& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}
    
    NGramTree& NGramTree::operator=(const NGramTree& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }


    void NGramTree::apply(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();

      pimpl->ngram_tree_score(state, states, edge, features);

      if (final)
	pimpl->ngram_tree_final_score(state, edge, features);
    }

    void NGramTree::apply_coarse(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {
    }
    void NGramTree::apply_predict(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {}
    void NGramTree::apply_scan(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       const int dot,
			       feature_set_type& features,
			       feature_set_type& estimates,
			       const bool final) const
    {}
    void NGramTree::apply_complete(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }


    void NGramTree::initialize()
    {
      pimpl->clear();
    }
    
    void NGramTree::assign(const size_type& id,
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
