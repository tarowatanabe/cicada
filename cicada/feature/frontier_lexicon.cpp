//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "frontier_lexicon.hpp"
#include "feature_builder.hpp"

#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"
#include "cicada/feature_vector_linear.hpp"

#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/array_power2.hpp"
#include "utils/small_vector.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>

namespace cicada
{
  namespace feature
  {


    class FrontierLexiconImpl : public utils::hashmurmur3<size_t>
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Lattice       lattice_type;
      typedef cicada::HyperGraph    hypergraph_type;

      typedef cicada::Cluster  cluster_type;
      typedef cicada::Stemmer  stemmer_type;
      
      typedef cicada::ClusterStemmer normalizer_type;
      typedef std::vector<normalizer_type, std::allocator<normalizer_type> > normalizer_set_type;

      typedef utils::hashmurmur3<size_t> hasher_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;

      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;

      typedef FeatureVectorLinear<feature_set_type::mapped_type> feature_linear_set_type;

      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;
      
      typedef symbol_type word_type;
            
      struct CacheNormalize
      {
	typedef utils::small_vector<word_type, std::allocator<word_type> > word_set_type;
	
	word_type::id_type word;
	word_set_type      normalized;
	
	CacheNormalize() : word(word_type::id_type(-1)), normalized() {}
      };
      typedef CacheNormalize cache_normalize_type;
      typedef utils::array_power2<cache_normalize_type, 1024 * 4, std::allocator<cache_normalize_type> > cache_normalize_set_type;
      
      struct cache_phrase_type
      {
	std::string frontier;
	phrase_type phrase;
	
	cache_phrase_type() : frontier(), phrase() {}
      };
      typedef utils::array_power2<cache_phrase_type, 1024 * 4, std::allocator<cache_phrase_type> > cache_phrase_set_type;
      
      struct cache_feature_type
      {
	word_type source;
	word_type target;
	feature_linear_set_type features;
	
	cache_feature_type() : source(), target(), features() {}
      };
      typedef utils::array_power2<cache_feature_type, 1024 * 4, std::allocator<cache_feature_type> > cache_feature_set_type;
      
      typedef FeatureBuilder feature_builder_type;

      FrontierLexiconImpl()
	: skip_sgml_tag(false), prefix("frontier-lexicon"), forced_feature(false),
	  attr_frontier_source("frontier-source"),
	  attr_frontier_target("frontier-target")
      { }
      
      struct skipper_epsilon
      {
	bool operator()(const symbol_type& word) const
	{
	  return word.is_non_terminal() || word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS;
	}
      };

      struct skipper_sgml
      {
	bool operator()(const symbol_type& word) const
	{
	  return word.is_non_terminal() || word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS || word.is_sgml_tag();
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
      
      struct __attribute_string : public boost::static_visitor<cicada::AttributeVector::string_type>
      {
	typedef cicada::AttributeVector attribute_set_type;
	
	static 
	const std::string& empty()
	{
	  static std::string __empty;
	  return __empty;
	}
	
	const attribute_set_type::string_type& operator()(const attribute_set_type::int_type& x) const { return empty(); }
	const attribute_set_type::string_type& operator()(const attribute_set_type::float_type& x) const { return empty(); }
	const attribute_set_type::string_type& operator()(const attribute_set_type::string_type& x) const { return x; }
      };
      
      template <typename Skipper>
      void lexicon_score(const edge_type& edge,
			 feature_set_type& features,
			 Skipper skipper)
      {
	attribute_set_type::const_iterator siter = edge.attributes.find(attr_frontier_source);
	attribute_set_type::const_iterator titer = edge.attributes.find(attr_frontier_target);
	
	if (siter == edge.attributes.end() || titer == edge.attributes.end()) return;
	
	const std::string& frontier_source = boost::apply_visitor(__attribute_string(), siter->second);
	const std::string& frontier_target = boost::apply_visitor(__attribute_string(), titer->second);
	
	if (frontier_source.empty() || frontier_target.empty()) return;
	
	const phrase_type& phrase_source = phrase(frontier_source, cache_phrase_source);
	const phrase_type& phrase_target = phrase(frontier_target, cache_phrase_target);
	
	if (phrase_source.empty() || phrase_target.empty()) return;

	phrase_type::const_iterator siter_begin = phrase_source.begin();
	phrase_type::const_iterator siter_end   = phrase_source.end();
	phrase_type::const_iterator titer_begin = phrase_target.begin();
	phrase_type::const_iterator titer_end   = phrase_target.end();
		
	for (phrase_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	  if (! skipper(*siter))
	    for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) 
	      if (! skipper(*titer)) {
		const size_t cache_pos = hasher_type::operator()(siter->id(), titer->id()) & (cache_features.size() - 1);
		cache_feature_type& cache = cache_features[cache_pos];
		
		if (cache.source != *siter || cache.target != *titer) {
		  cache.source = *siter;
		  cache.target = *titer;
		  
		  cache.features.clear();
		  
		  apply(*siter, *titer, cache.features);
		}
		
		features += cache.features;
	      }
      }
      
      template <typename Features>
      void apply(const word_type& source, const word_type& target, Features& features)
      {
	feature_builder.clear();
	feature_builder << prefix << ":" << source << "_" << target;
	
	if (forced_feature || feature_builder.exists())
	  features[feature_builder] += 1.0;
	
	if (! normalizers_source.empty()) {
	  const cache_normalize_type::word_set_type& normalized_source = normalize(source,  normalizers_source, cache_source);
	  
	  cache_normalize_type::word_set_type::const_iterator siter_end = normalized_source.end();
	  for (cache_normalize_type::word_set_type::const_iterator siter = normalized_source.begin(); siter != siter_end; ++ siter) {
	    
	    feature_builder.clear();
	    feature_builder << prefix << ":" << *siter << "_" << target;
	    
	    if (forced_feature || feature_builder.exists())
	      features[feature_builder] += 1.0;
	    
	    if (! normalizers_target.empty()) {
	      const cache_normalize_type::word_set_type& normalized_target = normalize(target, normalizers_target, cache_target);
	      
	      cache_normalize_type::word_set_type::const_iterator titer_end = normalized_target.end();
	      for (cache_normalize_type::word_set_type::const_iterator titer = normalized_target.begin(); titer != titer_end; ++ titer) {
		
		feature_builder.clear();
		feature_builder << prefix << ":" << *siter << "_" << *titer;
		
		if (forced_feature || feature_builder.exists())
		  features[feature_builder] += 1.0;
	      }
	    }
	  }
	}
	
	if (! normalizers_target.empty()) {
	  const cache_normalize_type::word_set_type& normalized_target = normalize(target, normalizers_target, cache_target);
	  
	  cache_normalize_type::word_set_type::const_iterator titer_end = normalized_target.end();
	  for (cache_normalize_type::word_set_type::const_iterator titer = normalized_target.begin(); titer != titer_end; ++ titer) {
	    feature_builder.clear();
	    feature_builder << prefix << ":" << source << "_" << *titer;
	    
	    if (forced_feature || feature_builder.exists())
	      features[feature_builder] += 1.0;
	  }
	}
      }

      const phrase_type& phrase(const std::string& frontier,
				cache_phrase_set_type& caches)
      {
	cache_phrase_type& cache = caches[hasher_type::operator()(frontier.begin(), frontier.end(), 0) & (caches.size() - 1)];
	
	if (cache.frontier != frontier) {
	  typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
	  
	  utils::piece frontier_piece(frontier);
	  tokenizer_type tokenizer(frontier_piece);
	  cache.phrase.assign(tokenizer.begin(), tokenizer.end());
	  cache.frontier = frontier;
	}
	
	return cache.phrase;
      }
      
      const cache_normalize_type::word_set_type& normalize(const word_type& word,
							   normalizer_set_type& normalizers,
							   cache_normalize_set_type& caches)
      {
	cache_normalize_type& cache = caches[hasher_type::operator()(word.id()) & (caches.size() - 1)];
	if (cache.word != word.id()) {
	  cache.word = word.id();
	  cache.normalized.clear();
	  
	  for (size_t i = 0; i != normalizers.size(); ++ i) {
	    const word_type normalized = normalizers[i](word);
	    if (word != normalized)
	      cache.normalized.push_back(normalized);
	  }
	}
	
	return cache.normalized;
      }
      
      void clear()
      {
      }
      
      void clear_cache()
      {
	cache_source.clear();
	cache_target.clear();

	cache_phrase_source.clear();
	cache_phrase_target.clear();

	cache_features.clear();
      }

      normalizer_set_type normalizers_source;
      normalizer_set_type normalizers_target;
      
      cache_normalize_set_type cache_source;
      cache_normalize_set_type cache_target;

      cache_phrase_set_type cache_phrase_source;
      cache_phrase_set_type cache_phrase_target;

      cache_feature_set_type cache_features;

      feature_builder_type feature_builder;
      
      bool skip_sgml_tag;
      
      std::string prefix;
      bool forced_feature;

      attribute_type attr_frontier_source;
      attribute_type attr_frontier_target;
    };
    
    FrontierLexicon::FrontierLexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "frontier-lexicon")
	throw std::runtime_error("this is not frontier lexicon feature: " + parameter);

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
	  std::cerr << "WARNING: unsupported parameter for frontier lexicon: " << piter->first << "=" << piter->second << std::endl;
      }
            
      std::auto_ptr<impl_type> lexicon_impl(new impl_type());

      lexicon_impl->normalizers_source.swap(normalizers_source);
      lexicon_impl->normalizers_target.swap(normalizers_target);
      
      lexicon_impl->skip_sgml_tag = skip_sgml_tag;
      lexicon_impl->prefix = (name.empty() ? std::string("frontier-lexicon") : name);

      base_type::__state_size = 0;
      base_type::__feature_name = (name.empty() ? std::string("frontier-lexicon") : name);
      base_type::__sparse_feature = true;
      
      pimpl = lexicon_impl.release();
    }
    
    FrontierLexicon::~FrontierLexicon() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    FrontierLexicon::FrontierLexicon(const FrontierLexicon& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {
      pimpl->clear_cache();
    }

    FrontierLexicon& FrontierLexicon::operator=(const FrontierLexicon& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      pimpl->clear_cache();
            
      return *this;
    }
    
    void FrontierLexicon::apply(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				const bool final) const
    {
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      feature_set_type feats;
      
      pimpl->lexicon_score(edge, feats);

      features.update(feats, static_cast<const std::string&>(base_type::feature_name()));
    }

    void FrontierLexicon::apply_coarse(state_ptr_type& state,
				       const state_ptr_set_type& states,
				       const edge_type& edge,
				       feature_set_type& features,
				       const bool final) const
    {
      apply(state, states, edge, features, final);
    }
    
    void FrontierLexicon::apply_predict(state_ptr_type& state,
					const state_ptr_set_type& states,
					const edge_type& edge,
					feature_set_type& features,
					const bool final) const
    {
      apply(state, states, edge, features, final);
    }

    
    void FrontierLexicon::apply_scan(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     const int dot,
				     feature_set_type& features,
				     const bool final) const
    {}
    void FrontierLexicon::apply_complete(state_ptr_type& state,
					 const state_ptr_set_type& states,
					 const edge_type& edge,
					 feature_set_type& features,
					 const bool final) const
    {}
    
  };
};
