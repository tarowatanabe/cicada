
#include <utility>
#include <memory>

#include "cicada/feature/ngram_tree.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"

#include "utils/indexed_set.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/lexical_cast.hpp"

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
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef rule_type::symbol_set_type phrase_type;
      
      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;
      
      struct node_pair_type
      {
	std::string node;
	std::string cluster;
	std::string prefix;
	std::string suffix;
	std::string digits;
	
	node_pair_type() : node(), cluster(), prefix(), suffix(), digits() {}
      };
      
      typedef utils::compact_trie_dense<symbol_type, node_pair_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
					std::allocator<std::pair<const symbol_type, node_pair_type> > > tree_map_type;

      typedef tree_map_type::id_type id_type;
      
      
      NGramTreeImpl()
	: cluster(0), stemmer_prefix(0), stemmer_suffix(0), stemmer_digits(0),
	  tree_map(symbol_type()),
	  forced_feature(false) {}
      
      virtual ~NGramTreeImpl() {}
      
      virtual void ngram_tree_score(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features) const = 0;
      virtual void ngram_tree_final_score(const state_ptr_type& state,
					  const edge_type& edge,
					  feature_set_type& features) const = 0;
      
      void clear()
      {
	tree_map.clear();
      }

      cluster_type* cluster;
      stemmer_type* stemmer_prefix;
      stemmer_type* stemmer_suffix;
      stemmer_type* stemmer_digits;
      
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
	
	const rule_type::symbol_set_type& phrase = extract_phrase(edge);
	
	if (states.empty()) {
	  // we do not add feature here, since we know nothing abount surrounding context...
	  symbol_type prefix = vocab_type::EPSILON;
	  symbol_type suffix = vocab_type::EPSILON;
	  
	  compute_bound(phrase.begin(), phrase.end(), prefix, suffix);
	  
	  id_type* context = reinterpret_cast<id_type*>(state);
	  //context[0] = tree_id(edge.rule->lhs, tree_id(prefix, tree_map.root()));
	  //context[1] = tree_id(edge.rule->lhs, tree_id(suffix, tree_map.root()));
	  
	  context[0] = tree_id(prefix, tree_map.root());
	  context[1] = tree_id(suffix, tree_map.root());
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
	  
	  phrase_span_set_type::const_iterator siter_begin = phrase_spans.begin();
	  phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	  for (phrase_span_set_type::const_iterator siter = siter_begin + 1; siter != siter_end; ++ siter) {
	    const phrase_span_type& span = *siter;
	    
	    // incase, we are working with non-synchronous parsing!
	    int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      antecedent_index = siter - (siter_begin + 1);
	    
	    const id_type* antecedent_context = reinterpret_cast<const id_type*>(states[antecedent_index]);
	    
	    //const id_type prefix_antecedent_id = antecedent_context[0];
	    //const id_type suffix_antecedent_id = antecedent_context[1];
	    const id_type prefix_antecedent_id = tree_id(*(span.first - 1), antecedent_context[0]);
	    const id_type suffix_antecedent_id = tree_id(*(span.first - 1), antecedent_context[1]);
	    
	    symbol_type prefix_next;
	    symbol_type suffix_next;
	    compute_bound(span.first, span.second, prefix_next, suffix_next);
	    
	    const id_type prefix_next_id = (prefix_next.empty() ? tree_map.root() : tree_id(prefix_next, tree_map.root()));
	    const id_type suffix_next_id = (suffix_next.empty() ? tree_map.root() : tree_id(suffix_next, tree_map.root()));
	    
	    if (! tree_map.is_root(suffix_id))
	      apply_feature(features, edge.rule->lhs, suffix_id, prefix_antecedent_id);
	    
	    if (tree_map.is_root(prefix_id))
	      prefix_id = prefix_antecedent_id;
	    
	    if (! tree_map.is_root(prefix_next_id)) {
	      apply_feature(features, edge.rule->lhs, suffix_antecedent_id, prefix_next_id);
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
	  
	  //context[0] = tree_id(edge.rule->lhs, prefix_id);
	  //context[1] = tree_id(edge.rule->lhs, suffix_id);
	  context[0] = prefix_id;
	  context[1] = suffix_id;
	}
      }

      virtual void ngram_tree_final_score(const state_ptr_type& state,
					  const edge_type& edge,
					  feature_set_type& features) const
      {
	const id_type* antecedent_context = reinterpret_cast<const id_type*>(state);
	
	const id_type prefix_antecedent_id = antecedent_context[0];
	const id_type suffix_antecedent_id = antecedent_context[1];
	
	apply_feature(features, edge.rule->lhs, tree_id(vocab_type::BOS, tree_map.root()), prefix_antecedent_id);
	apply_feature(features, edge.rule->lhs, suffix_antecedent_id, tree_id(vocab_type::EOS, tree_map.root()));
      }

      
      id_type tree_id(const symbol_type& node, const id_type parent) const
      {
	tree_map_type& __tree_map = const_cast<tree_map_type&>(tree_map);

	const id_type id = __tree_map.insert(parent, node);
	
	if (__tree_map[id].node.empty()) {
	  if (__tree_map.is_root(parent)) {
	    __tree_map[id].node = node;
	    if (cluster)
	      __tree_map[id].cluster = cluster->operator[](node);
	    if (stemmer_prefix)
	      __tree_map[id].prefix = stemmer_prefix->operator[](node);
	    if (stemmer_suffix)
	      __tree_map[id].suffix = stemmer_suffix->operator[](node);
	    if (stemmer_digits)
	      __tree_map[id].digits = stemmer_digits->operator[](node);
	  } else {
	    __tree_map[id].node = compose_path(node, __tree_map[parent].node);
	    if (cluster)
	      __tree_map[id].cluster = compose_path(node, __tree_map[parent].cluster);
	    if (stemmer_prefix)
	      __tree_map[id].prefix = compose_path(node, __tree_map[parent].prefix);
	    if (stemmer_suffix)
	      __tree_map[id].suffix = compose_path(node, __tree_map[parent].suffix);
	    if (stemmer_digits)
	      __tree_map[id].digits = compose_path(node, __tree_map[parent].digits);
	  }
	}
	return id;
      }

      void apply_feature(feature_set_type& features, const std::string& node, const id_type& prev, const id_type& next) const
      {
	const node_pair_type& prev_node = tree_map[prev];
	const node_pair_type& next_node = tree_map[next];
	
	if (cluster && (prev_node.node != prev_node.cluster || next_node.node != next_node.cluster)) {
	  const std::string name = feature_name(node, prev_node.cluster, next_node.cluster);
	  if (forced_feature || feature_set_type::feature_type::exists(name))
	    features[name] += 1.0;
	}
	
	if (stemmer_prefix && (prev_node.node != prev_node.prefix || next_node.node != next_node.prefix)) {
	  const std::string name = feature_name(node, prev_node.prefix, next_node.prefix);
	  if (forced_feature || feature_set_type::feature_type::exists(name))
	    features[name] += 1.0;
	}

	if (stemmer_suffix && (prev_node.node != prev_node.suffix || next_node.node != next_node.suffix)) {
	  const std::string name = feature_name(node, prev_node.suffix, next_node.suffix);
	  if (forced_feature || feature_set_type::feature_type::exists(name))
	    features[name] += 1.0;
	}
	
	if (stemmer_digits && (prev_node.node != prev_node.digits || next_node.node != next_node.digits)) {
	  const std::string name = feature_name(node, prev_node.digits, next_node.digits);
	  if (forced_feature || feature_set_type::feature_type::exists(name))
	    features[name] += 1.0;
	}
	
	
	const std::string name = feature_name(node, prev_node.node, next_node.node);
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;
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
	return Extract::feature_prefix +  compose_tree(node, prev, next);
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
	  std::cerr << "WARNING: unsupported parameter for ngram-tree: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (source && target)
	throw std::runtime_error("both source and target?");
      if (! source && ! target)
	throw std::runtime_error("what side are you going to use?");

      if (stemmer_prefix_size < 0)
	throw std::runtime_error("negative prefix size?");
      if (stemmer_suffix_size < 0)
	throw std::runtime_error("negative suffix size?");
      
      std::auto_ptr<impl_type> ngram_tree_impl(source
					       ? dynamic_cast<impl_type*>(new __NGramTreeImpl<__ngram_tree_extract_source>())
					       : dynamic_cast<impl_type*>(new __NGramTreeImpl<__ngram_tree_extract_target>()));

      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file: " + cluster_path.file_string());
	
	ngram_tree_impl->cluster = &cicada::Cluster::create(cluster_path);
      }
      
      if (stemmer_prefix_size > 0)
	ngram_tree_impl->stemmer_prefix = &cicada::Stemmer::create("prefix:size=" + boost::lexical_cast<std::string>(stemmer_prefix_size));
      
      if (stemmer_suffix_size > 0)
	ngram_tree_impl->stemmer_suffix = &cicada::Stemmer::create("suffix:size=" + boost::lexical_cast<std::string>(stemmer_suffix_size));

      if (stemmer_digits)
	ngram_tree_impl->stemmer_digits = &cicada::Stemmer::create("digits");

      
      // non-terminal + two neighbouring symbols + span-size
      base_type::__state_size = sizeof(impl_type::id_type) * 2;
      base_type::__feature_name = std::string("ngram-tree-") + (source ? "source" : "target");
      base_type::__sparse_feature = true;
      
      pimpl = ngram_tree_impl.release();
    }
    
    NGramTree::~NGramTree() { std::auto_ptr<impl_type> tmp(pimpl); }

    
    NGramTree::NGramTree(const NGramTree& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(0)
    {
      typedef __NGramTreeImpl<__ngram_tree_extract_source> ngram_tree_source_type;
      typedef __NGramTreeImpl<__ngram_tree_extract_target> ngram_tree_target_type;
      
      if (dynamic_cast<const ngram_tree_source_type*>(x.pimpl))
	pimpl = new ngram_tree_source_type(*dynamic_cast<const ngram_tree_source_type*>(x.pimpl));
      else
	pimpl = new ngram_tree_target_type(*dynamic_cast<const ngram_tree_target_type*>(x.pimpl));
    }
    
    NGramTree& NGramTree::operator=(const NGramTree& x)
    {
      typedef __NGramTreeImpl<__ngram_tree_extract_source> ngram_tree_source_type;
      typedef __NGramTreeImpl<__ngram_tree_extract_target> ngram_tree_target_type;

      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      std::auto_ptr<impl_type> tmp(pimpl);
      
      if (dynamic_cast<const ngram_tree_source_type*>(x.pimpl))
	pimpl = new ngram_tree_source_type(*dynamic_cast<const ngram_tree_source_type*>(x.pimpl));
      else
	pimpl = new ngram_tree_target_type(*dynamic_cast<const ngram_tree_target_type*>(x.pimpl));

      return *this;
    }


    template <typename FeaturePrefix, typename Feature>
    inline
    bool equal_prefix(const FeaturePrefix& prefix, const Feature& x)
    {
      return x.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), x.begin());
    }
    
    void NGramTree::apply(state_ptr_type& state,
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
  };
};
