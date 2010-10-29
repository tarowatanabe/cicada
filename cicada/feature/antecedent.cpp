
#include <utility>
#include <memory>

#include "cicada/feature/antecedent.hpp"
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
    
    class AntecedentImpl
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
      
      typedef utils::compact_trie_dense<symbol_type, std::string, boost::hash<symbol_type>, std::equal_to<symbol_type>,
					std::allocator<std::pair<const symbol_type, std::string> > > tree_map_type;
      
      typedef tree_map_type::id_type id_type;
      
      AntecedentImpl()
	: cluster(0), stemmer_prefix(0), stemmer_suffix(0), stemmer_digits(0),
	  tree_map(symbol_type()),
	  forced_feature(false) {}
      
      
      void clear()
      {
	tree_map.clear();
      }

      cluster_type* cluster;
      stemmer_type* stemmer_prefix;
      stemmer_type* stemmer_suffix;
      stemmer_type* stemmer_digits;
      
      tree_map_type  tree_map;

      bool forced_feature;
      
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

	  phrase_type::const_iterator piter_end = phrase.end();
	  for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter)
	    if (*piter != vocab_type::EPSILON) {
	      if (prefix == vocab_type::EPSILON)
		prefix = *piter;
	      suffix = *piter;
	      ++ span_size;
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
	if (cluster) {
	  const symbol_type prefix_cluster = cluster->operator[](prefix);
	  const symbol_type suffix_cluster = cluster->operator[](suffix);

	  if (prefix_cluster != prefix || suffix_cluster != suffix) {
	    const std::string name = feature_name(node, antecedent, prefix_cluster, suffix_cluster, span_size);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}
	
	if (stemmer_prefix) {
	  const symbol_type prefix_stemmed = stemmer_prefix->operator[](prefix);
	  const symbol_type suffix_stemmed = stemmer_prefix->operator[](suffix);

	  if (prefix_stemmed != prefix || suffix_stemmed != suffix) {
	    const std::string name = feature_name(node, antecedent, prefix_stemmed, suffix_stemmed, span_size);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}

	if (stemmer_suffix) {
	  const symbol_type prefix_stemmed = stemmer_suffix->operator[](prefix);
	  const symbol_type suffix_stemmed = stemmer_suffix->operator[](suffix);

	  if (prefix_stemmed != prefix || suffix_stemmed != suffix) {
	    const std::string name = feature_name(node, antecedent, prefix_stemmed, suffix_stemmed, span_size);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}

	if (stemmer_digits) {
	  const symbol_type prefix_stemmed = stemmer_digits->operator[](prefix);
	  const symbol_type suffix_stemmed = stemmer_digits->operator[](suffix);

	  if (prefix_stemmed != prefix || suffix_stemmed != suffix) {
	    const std::string name = feature_name(node, antecedent, prefix_stemmed, suffix_stemmed, span_size);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}

	
	const std::string name = feature_name(node, antecedent, prefix, suffix, span_size);
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;
      }
      
      const std::string feature_name(const std::string& node,
				     const std::string& antecedent,
				     const std::string& prefix,
				     const std::string& suffix,
				     const int span_size) const
      {
	return "antecedent:" + node + antecedent + '|' + prefix + '|' + suffix + '|' + boost::lexical_cast<std::string>(span_size);
      }
    };
    
    
    Antecedent::Antecedent(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "antecedent")
	throw std::runtime_error("is this really antecedent feature function? " + parameter);
      
      int stemmer_prefix_size = 0;
      int stemmer_suffix_size = 0;
      bool stemmer_digits = false;
      
      boost::filesystem::path cluster_path;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "cluster") == 0)
	  cluster_path = piter->second;
	else if (strcasecmp(piter->first.c_str(), "prefix") == 0)
	  stemmer_prefix_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "suffix") == 0)
	  stemmer_suffix_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "digits") == 0)
	  stemmer_digits = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for antecedent: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (stemmer_prefix_size < 0)
	throw std::runtime_error("negative prefix size?");
      if (stemmer_suffix_size < 0)
	throw std::runtime_error("negative suffix size?");

      std::auto_ptr<impl_type> antecedent_impl(new impl_type());

      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file: " + cluster_path.file_string());
	
	antecedent_impl->cluster = &cicada::Cluster::create(cluster_path);
      }
      
      if (stemmer_prefix_size > 0)
	antecedent_impl->stemmer_prefix = &cicada::Stemmer::create("prefix:size=" + boost::lexical_cast<std::string>(stemmer_prefix_size));
      
      if (stemmer_suffix_size > 0)
	antecedent_impl->stemmer_suffix = &cicada::Stemmer::create("suffix:size=" + boost::lexical_cast<std::string>(stemmer_suffix_size));

      if (stemmer_digits)
	antecedent_impl->stemmer_digits = &cicada::Stemmer::create("digits");
      
      // antecedent conext + terminal-boundary + span-size
      base_type::__state_size = sizeof(impl_type::id_type) + sizeof(symbol_type) * 2 + sizeof(int);
      base_type::__feature_name = std::string("antecedent");
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
  };
};
