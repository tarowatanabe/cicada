
#include <utility>
#include <memory>

#include "cicada/feature/ngram_tree.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"

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
      typedef cicada::Cluster  cluster_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef rule_type::symbol_set_type phrase_type;
      
      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;

      typedef std::pair<std::string, std::string> node_pair_type;
      
      typedef utils::compact_trie<symbol_type, node_pair_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				  std::allocator<std::pair<const symbol_type, node_pair_type> > > tree_map_type;

      typedef tree_map_type::id_type id_type;
      
      
      NGramTreeImpl() : forced_feature(false) {}
      NGramTreeImpl(const NGramTreeImpl& x) : cluster(x.cluster), forced_feature(x.forced_feature) {}
      NGramTreeImpl& operator=(const NGramTreeImpl& x)
      {
	cluster = x.cluster;
	forced_feature = x.forced_feature;
	return *this;
      }

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

      cluster_type cluster;
      
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
	  symbol_type prefix = vocab_type::EPSILON;
	  symbol_type suffix = vocab_type::EPSILON;
	  
	  compute_bound(phrase.begin(), phrase.end(), prefix, suffix);
	  
	  id_type* context = reinterpret_cast<id_type*>(state);
	  context[0] = tree_id(edge.rule->lhs, tree_id(prefix, tree_map.root()));
	  context[1] = tree_id(edge.rule->lhs, tree_id(suffix, tree_map.root()));
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

	    const id_type prefix_antecedent_id = antecedent_context[0];
	    const id_type suffix_antecedent_id = antecedent_context[1];
	    
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
	  
	  context[0] = tree_id(edge.rule->lhs, prefix_id);
	  context[1] = tree_id(edge.rule->lhs, suffix_id);
	}
      }

      virtual void ngram_tree_final_score(const state_ptr_type& state,
					  feature_set_type& features) const
      {
	const id_type* antecedent_context = reinterpret_cast<const id_type*>(state);
	
	const id_type prefix_antecedent_id = antecedent_context[0];
	const id_type suffix_antecedent_id = antecedent_context[1];
	
	apply_feature(features, vocab_type::GOAL, tree_id(vocab_type::BOS, tree_map.root()), prefix_antecedent_id);
	apply_feature(features, vocab_type::GOAL, tree_id(vocab_type::EOS, tree_map.root()), suffix_antecedent_id);
      }

      
      id_type tree_id(const symbol_type& node, const id_type parent) const
      {
	tree_map_type& __tree_map = const_cast<tree_map_type&>(tree_map);

	const id_type id = __tree_map.insert(parent, node);
	
	if (__tree_map[id].first.empty()) {
	  if (__tree_map.is_root(parent)) {
	    __tree_map[id].first  = node;
	    __tree_map[id].second = cluster[node];
	  } else {
	    __tree_map[id].first  = compose_path(node, __tree_map[parent].first);
	    __tree_map[id].second = compose_path(node, __tree_map[parent].second);
	  }
	}
	return id;
      }

      void apply_feature(feature_set_type& features, const std::string& node, const id_type& prev, const id_type& next) const
      {
	const node_pair_type& prev_node = tree_map[prev];
	const node_pair_type& next_node = tree_map[next];
	
	if (! cluster.empty() && prev_node.first != prev_node.second && next_node.first != next_node.second) {
	  const std::string name = feature_name(node, prev_node.second, next_node.second);
	  if (forced_feature || feature_set_type::feature_type::exists(name))
	    features[name] += 1.0;
	}
	
	const std::string name = feature_name(node, prev_node.first, next_node.first);
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
	return node + "(" + prev + ")(" + next + ')';
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
	else
	  std::cerr << "WARNING: unsupported parameter for ngram-tree: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (source && target)
	throw std::runtime_error("both source and target?");
      if (! source && ! target)
	throw std::runtime_error("what side are you going to use?");
      
      std::auto_ptr<impl_type> ngram_tree_impl(source
					       ? dynamic_cast<impl_type*>(new __NGramTreeImpl<__ngram_tree_extract_source>())
					       : dynamic_cast<impl_type*>(new __NGramTreeImpl<__ngram_tree_extract_target>()));

      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file" + cluster_path.file_string());
	
	ngram_tree_impl->cluster = cicada::Cluster(cluster_path);
      }

      
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
	pimpl = new ngram_tree_source_type();
      else
	pimpl = new ngram_tree_target_type();
    }
    
    NGramTree& NGramTree::operator=(const NGramTree& x)
    {
      typedef __NGramTreeImpl<__ngram_tree_extract_source> ngram_tree_source_type;
      typedef __NGramTreeImpl<__ngram_tree_extract_target> ngram_tree_target_type;

      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      std::auto_ptr<impl_type> tmp(pimpl);
      
      if (dynamic_cast<const ngram_tree_source_type*>(x.pimpl))
	pimpl = new ngram_tree_source_type();
      else
	pimpl = new ngram_tree_target_type();
      
      return *this;
    }


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
