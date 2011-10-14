//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "sparse_lexicon.hpp"

#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"

#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/alloc_vector.hpp"

#include <google/dense_hash_set>
#include <google/dense_hash_map>

namespace cicada
{
  namespace feature
  {

    class SparseLexiconImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      typedef cicada::Lattice       lattice_type;
      typedef cicada::HyperGraph    hypergraph_type;

      typedef cicada::Cluster  cluster_type;
      typedef cicada::Stemmer  stemmer_type;
      
      typedef cicada::ClusterStemmer normalizer_type;
      typedef std::vector<normalizer_type, std::allocator<normalizer_type> > normalizer_set_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;

      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;

      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;
      
      typedef symbol_type word_type;
      typedef std::pair<word_type, word_type> word_pair_type;
      
      typedef google::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type> > word_set_type;

      typedef utils::alloc_vector<feature_set_type, std::allocator<feature_set_type> > cache_set_type;
      
      SparseLexiconImpl()
	: uniques(), words(), caches(), skip_sgml_tag(false), prefix("sparse-lexicon"), forced_feature(false)
      { uniques.set_empty_key(word_type());  }

      struct skipper_epsilon
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS;
	}
      };

      struct skipper_sgml
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS || word.is_sgml_tag();
	}
      };
      
      void lexicon_score(const edge_type& edge,
			 feature_set_type& features)
      {
	if (skip_sgml_tag)
	  lexicon_score(edge, features, skipper_sgml());
	else
	  lexicon_score(edge, features, skipper_epsilon());
      }
      
      
      template <typename Skipper>
      void lexicon_score(const edge_type& edge,
			 feature_set_type& features,
			 Skipper skipper)
      {
	const phrase_type& phrase = edge.rule->rhs;
	
	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter) 
	  if (piter->is_terminal() && ! skipper(*piter)) {

	    if (! caches.exists(piter->id())) {
	      feature_set_type& features = caches[piter->id()];
	      features.clear();
	      
	      sentence_type::const_iterator witer_end = words.end();
	      for (sentence_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer) {
		const std::string name = prefix + ":" + static_cast<const std::string&>(*witer) + "_" + static_cast<const std::string&>(*piter);
		
		if (forced_feature || feature_type::exists(name))
		  features[name] += 1.0;
		
		for (size_t i = 0; i != normalizers_target.size(); ++ i) {
		  const word_type normalized = normalizers_target[i](*piter);
		  if (normalized != *piter) {
		    const std::string name = prefix + ":" + static_cast<const std::string&>(*witer) + "_" + static_cast<const std::string&>(normalized);
		    
		    if (forced_feature || feature_type::exists(name))
		      features[name] += 1.0;
		  }
		}
	      }
	    }
	    
	    features += caches[piter->id()];
	  }
      }

      void assign(const lattice_type& lattice)
      {
	if (skip_sgml_tag)
	  assign(lattice, skipper_sgml());
	else
	  assign(lattice, skipper_epsilon());
      }
      
      void assign(const hypergraph_type& graph)
      {
	if (skip_sgml_tag)
	  assign(graph, skipper_sgml());
	else
	  assign(graph, skipper_epsilon());
      }
      
      template <typename Skipper>
      void assign(const lattice_type& lattice,
		  Skipper skipper)
      {
	uniques.clear();
	caches.clear();
	
	lattice_type::const_iterator liter_end = lattice.end();
	for (lattice_type::const_iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	  lattice_type::arc_set_type::const_iterator aiter_end = liter->end();
	  for (lattice_type::arc_set_type::const_iterator aiter = liter->begin(); aiter != aiter_end; ++ aiter)
	    if (! skipper(aiter->label))
	      uniques.insert(aiter->label);
	}

	words.clear();
	word_set_type::const_iterator uiter_end = uniques.end();
	for (word_set_type::const_iterator uiter = uniques.begin(); uiter != uiter_end; ++ uiter) {
	  words.push_back(*uiter);
	  
	  for (size_t i = 0; i != normalizers_source.size(); ++ i) {
	    const word_type normalized = normalizers_source[i](*uiter);
	    if (normalized != *uiter)
	      words.push_back(normalized);
	  }
	}
      }
      
      template <typename Skipper>
      void assign(const hypergraph_type& forest,
		  Skipper skipper)
      {
	uniques.clear();
	caches.clear();
	
	hypergraph_type::edge_set_type::const_iterator eiter_end = forest.edges.end();
	for (hypergraph_type::edge_set_type::const_iterator eiter = forest.edges.begin(); eiter != eiter_end; ++ eiter)
	  if (eiter->rule) {
	    const rule_type& rule = *(eiter->rule);
	    
	    rule_type::symbol_set_type::const_iterator siter_end = rule.rhs.end();
	    for (rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter)
	      if (siter->is_terminal() && ! skipper(*siter))
		uniques.insert(*siter);
	  }
	
	words.clear();
	word_set_type::const_iterator uiter_end = uniques.end();
	for (word_set_type::const_iterator uiter = uniques.begin(); uiter != uiter_end; ++ uiter) {
	  words.push_back(*uiter);
	  
	  for (size_t i = 0; i != normalizers_source.size(); ++ i) {
	    const word_type normalized = normalizers_source[i](*uiter);
	    if (normalized != *uiter)
	      words.push_back(normalized);
	  }
	}
      }
      
      void clear()
      {
	words.clear();
	caches.clear();
      }

      normalizer_set_type normalizers_source;
      normalizer_set_type normalizers_target;
      
      word_set_type  uniques;
      sentence_type  words;
      cache_set_type caches;
      
      bool skip_sgml_tag;
      
      std::string prefix;
      bool forced_feature;
    };
    
    SparseLexicon::SparseLexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "sparse-lexicon")
	throw std::runtime_error("this is not sparse lexicon feature: " + parameter);
      
      impl_type::normalizer_set_type normalizers_source;
      impl_type::normalizer_set_type normalizers_target;
      
      bool skip_sgml_tag = false;
      
      std::string name;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "cluster-source") {
	  if (! boost::filesystem::exists(piter->second))
	    throw std::runtime_error("no cluster file: " + piter->second);
	  
	  normalizers_source.push_back(impl_type::normalizer_type(&cicada::Cluster::create(piter->second)));
	} else if (utils::ipiece(piter->first) == "cluster-target") {
	  if (! boost::filesystem::exists(piter->second))
	    throw std::runtime_error("no cluster file: " + piter->second);
	  
	  normalizers_target.push_back(impl_type::normalizer_type(&cicada::Cluster::create(piter->second)));
	} else if (utils::ipiece(piter->first) == "stemmer-source")
	  normalizers_source.push_back(impl_type::normalizer_type(&cicada::Stemmer::create(piter->second)));
	else if (utils::ipiece(piter->first) == "stemmer-target")
	  normalizers_target.push_back(impl_type::normalizer_type(&cicada::Stemmer::create(piter->second)));
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for sparse lexicon: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> lexicon_impl(new impl_type());

      lexicon_impl->normalizers_source.swap(normalizers_source);
      lexicon_impl->normalizers_target.swap(normalizers_target);
      
      lexicon_impl->skip_sgml_tag = skip_sgml_tag;
      lexicon_impl->prefix = (name.empty() ? std::string("sparse-lexicon") : name);
      
      // two-side context + length + counts-id 
      base_type::__state_size = 0;
      base_type::__feature_name = (name.empty() ? std::string("sparse-lexicon") : name);
      base_type::__sparse_feature = true;
      
      pimpl = lexicon_impl.release();
    }
    
    SparseLexicon::~SparseLexicon() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    SparseLexicon::SparseLexicon(const SparseLexicon& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    SparseLexicon& SparseLexicon::operator=(const SparseLexicon& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void SparseLexicon::apply(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));

      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->lexicon_score(edge, features);
    }

    void SparseLexicon::apply_coarse(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     feature_set_type& features,
				     feature_set_type& estimates,
				     const bool final) const
    {
      
    }
    
    void SparseLexicon::apply_predict(state_ptr_type& state,
				      const state_ptr_set_type& states,
				      const edge_type& edge,
				      feature_set_type& features,
				      feature_set_type& estimates,
				      const bool final) const
    {}
    void SparseLexicon::apply_scan(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   const int dot,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {}
    void SparseLexicon::apply_complete(state_ptr_type& state,
				       const state_ptr_set_type& states,
				       const edge_type& edge,
				       feature_set_type& features,
				       feature_set_type& estimates,
				       const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    
    void SparseLexicon::assign(const size_type& id,
			       const hypergraph_type& hypergraph,
			       const lattice_type& lattice,
			       const span_set_type& spans,
			       const sentence_set_type& targets,
			       const ngram_count_set_type& ngram_counts)
    {
      //
      // how do we assign lexion feature from hypergraph...???
      // we assume that the lattice is always filled with the source-word...!
      //
      
      pimpl->clear();
      
      if (! lattice.empty())
	pimpl->assign(lattice);
    }
    
  };
};
