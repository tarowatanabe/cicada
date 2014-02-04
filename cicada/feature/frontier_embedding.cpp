//
//  Copyright(C) 2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <memory>

#include "cicada/feature/frontier_embedding.hpp"
#include "cicada/word_embedding.hpp"
#include "cicada/parameter.hpp"
#include "cicada/symbol_vector.hpp"

#include "utils/array_power2.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"
#include "utils/random_seed.hpp"
#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>

namespace cicada
{
  namespace feature
  {
    class FrontierEmbeddingImpl : public utils::hashmurmur3<size_t>
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef utils::hashmurmur3<size_t> hasher_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Vocab  vocab_type;
      
      typedef cicada::WordEmbedding word_embedding_type;

      typedef word_embedding_type::parameter_type parameter_type;
      typedef word_embedding_type::tensor_type    tensor_type;
      typedef word_embedding_type::matrix_type    matrix_type;

      typedef boost::filesystem::path path_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;

      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;

      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
      typedef std::vector<feature_type, std::allocator<feature_type> > feature_name_set_type;
      
      typedef rule_type::word_type word_type;

      typedef rule_type::symbol_set_type phrase_type;
      
      struct cache_phrase_type
      {
	std::string frontier;
	phrase_type phrase;
	
	cache_phrase_type() : frontier(), phrase() {}
      };
      typedef utils::array_power2<cache_phrase_type, 1024 * 4, std::allocator<cache_phrase_type> > cache_phrase_set_type;

      struct __attribute_string : public boost::static_visitor<const cicada::AttributeVector::string_type&>
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
      
    public:
      FrontierEmbeddingImpl(const path_type& path_source,
			    const path_type& path_target,
			    const std::string& name)
	: embedding_source_(&word_embedding_type::create(path_source)),
	  embedding_target_(&word_embedding_type::create(path_target)),
	  path_source_(path_source),
	  path_target_(path_target),
	  attr_frontier_source_("frontier-source"),
	  attr_frontier_target_("frontier-target")
      {
	initialize(name);
      }
      
      FrontierEmbeddingImpl(const FrontierEmbeddingImpl& x)
	: embedding_source_(&word_embedding_type::create(x.path_source_)),
	  embedding_target_(&word_embedding_type::create(x.path_target_)),
	  path_source_(x.path_source_),
	  path_target_(x.path_target_),
	  skip_sgml_tag_(x.skip_sgml_tag_),
	  feature_names_source_(x.feature_names_source_),
	  feature_names_target_(x.feature_names_target_),
	  attr_frontier_source_("frontier-source"),
	  attr_frontier_target_("frontier-target")
      { }
      
      FrontierEmbeddingImpl& operator=(const FrontierEmbeddingImpl& x)
      {
	embedding_source_ = &word_embedding_type::create(x.path_source_);
	embedding_target_ = &word_embedding_type::create(x.path_target_);
	
	path_source_ = x.path_source_;
	path_target_ = x.path_target_;
	
	skip_sgml_tag_ = x.skip_sgml_tag_;
	
	feature_names_source_ = x.feature_names_source_;
	feature_names_target_ = x.feature_names_target_;
	
	cache_source_.clear();
	cache_target_.clear();
	
	return *this;
      }
      
      void initialize(const std::string& name)
      {
	if (! embedding_source_)
	  throw std::runtime_error("no source embedding?");

	if (! embedding_target_)
	  throw std::runtime_error("no target embedding?");
	
	feature_names_source_.clear();
	feature_names_target_.clear();
	
	for (size_type i = 0; i != embedding_source_->dimension(); ++ i)
	  feature_names_source_.push_back(name + ":source:" + utils::lexical_cast<std::string>(i));
	
	for (size_type i = 0; i != embedding_target_->dimension(); ++ i)
	  feature_names_target_.push_back(name + ":target:" + utils::lexical_cast<std::string>(i));
      }
      
      struct extract_word
      {
	const symbol_type& operator()(const symbol_type& word) const
	{
	  return word;
	}
      };

      struct skipper_epsilon
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON;
	}
      };
      
      struct skipper_sgml
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || (word != vocab_type::BOS && word != vocab_type::EOS && word.is_sgml_tag());
	}
      };
      
      void embedding_score(const edge_type& edge,
			   feature_set_type& features) const
      {
	if (skip_sgml_tag_)
	  embedding_score(edge, features, extract_word(), skipper_sgml());
	else
	  embedding_score(edge, features, extract_word(), skipper_epsilon());
      }

      
      template <typename Extract, typename Skipper>
      void embedding_score(const edge_type& edge,
			   feature_set_type& features,
			   Extract extract,
			   Skipper skipper) const
      {
	attribute_set_type::const_iterator siter = edge.attributes.find(attr_frontier_source_);
	attribute_set_type::const_iterator titer = edge.attributes.find(attr_frontier_target_);

	if (siter != edge.attributes.end() || titer != edge.attributes.end()) {
	  const std::string& frontier_source = (siter != edge.attributes.end()
						? boost::apply_visitor(__attribute_string(), siter->second)
						: frontier_tmp_);
	  const std::string& frontier_target = (titer != edge.attributes.end()
						? boost::apply_visitor(__attribute_string(), titer->second)
						: frontier_tmp_);
	  
	  const phrase_type& phrase_source = cache_phrase(frontier_source, cache_source_, skipper);
	  const phrase_type& phrase_target = cache_phrase(frontier_target, cache_target_, skipper);
	  
	  if (! phrase_source.empty()) {
	    tensor_type layer = tensor_type::Zero(embedding_source_->dimension(), 0);
	    
	    phrase_type::const_iterator piter_end = phrase_source.end();
	    for (phrase_type::const_iterator piter = phrase_source.begin(); piter != piter_end; ++ piter)
	      layer += embedding_source_->operator()(*piter);

	    for (size_type i = 0; i != embedding_source_->dimension(); ++ i)
	      if (layer(i, 0) != parameter_type(0))
		features[feature_names_source_[i]] = layer(i, 0);
	      else
		features.erase(feature_names_source_[i]);
	  }
	  
	  if (! phrase_target.empty()) {
	    tensor_type layer = tensor_type::Zero(embedding_target_->dimension(), 0);
	    
	    phrase_type::const_iterator piter_end = phrase_target.end();
	    for (phrase_type::const_iterator piter = phrase_target.begin(); piter != piter_end; ++ piter)
	      layer += embedding_target_->operator()(*piter);

	    for (size_type i = 0; i != embedding_target_->dimension(); ++ i)
	      if (layer(i, 0) != parameter_type(0))
		features[feature_names_target_[i]] = layer(i, 0);
	      else
		features.erase(feature_names_target_[i]);
	  }
	}
      }
      
      template <typename Skipper>
      const phrase_type& cache_phrase(const std::string& frontier,
				      const cache_phrase_set_type& caches,
				      Skipper skipper) const
      {
	if (frontier.empty()) return phrase_tmp_;
	
	const size_type cache_pos = hasher_type::operator()(frontier.begin(), frontier.end(), 0)& (caches.size() - 1);
	
	cache_phrase_type& cache = const_cast<cache_phrase_type&>(caches[cache_pos]);
	
	if (cache.frontier != frontier) {
	  typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
	  
	  utils::piece frontier_piece(frontier);
	  tokenizer_type tokenizer(frontier_piece);
	  
	  cache.frontier = frontier;
	  
	  cache.phrase.clear();
	  tokenizer_type::iterator titer_end = tokenizer.end();
	  for (tokenizer_type::iterator titer = tokenizer.begin(); titer != titer_end; ++ titer) {
	    const symbol_type word = *titer;
	    
	    if (! word.is_non_terminal() && ! skipper(word))
	      cache.phrase.push_back(word);
	  }
	}
	
	return cache.phrase;
      }

    public:
      const word_embedding_type* embedding_source_;
      const word_embedding_type* embedding_target_;
      
      path_type path_source_;
      path_type path_target_;
      
      bool skip_sgml_tag_;

      std::string   frontier_tmp_;
      phrase_type   phrase_tmp_;

      cache_phrase_set_type cache_source_;
      cache_phrase_set_type cache_target_;
      
      // names...
      feature_name_set_type feature_names_source_;
      feature_name_set_type feature_names_target_;

      attribute_type attr_frontier_source_;
      attribute_type attr_frontier_target_;
    };
    
    FrontierEmbedding::FrontierEmbedding(const std::string& parameter)
      : pimpl(0)
    {
      typedef boost::filesystem::path path_type;
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "frontier-embedding")
	throw std::runtime_error("is this really embedding feature function? " + parameter);
      
      path_type   path_source;
      path_type   path_target;
      bool        skip_sgml_tag = false;
      std::string name;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "source")
	  path_source = piter->second;
	else if (utils::ipiece(piter->first) == "target")
	  path_target = piter->second;
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for ngram: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (! boost::filesystem::exists(path_source))
	throw std::runtime_error("no source embedding file? " + path_source.string());

      if (! boost::filesystem::exists(path_target))
	throw std::runtime_error("no target embedding file? " + path_target.string());
      
      if (name.empty())
	name = "frontier-embedding";
      
      std::auto_ptr<impl_type> embedding_impl(new impl_type(path_source, path_target, name));
      
      embedding_impl->skip_sgml_tag_ = skip_sgml_tag;
      
      base_type::__state_size   = 0;
      base_type::__feature_name = name;
      
      pimpl = embedding_impl.release();
    }
    
    FrontierEmbedding::~FrontierEmbedding()
    {
      std::auto_ptr<impl_type> tmp(pimpl);
    }
    
    FrontierEmbedding::FrontierEmbedding(const FrontierEmbedding& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}
    
    FrontierEmbedding& FrontierEmbedding::operator=(const FrontierEmbedding& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    
    void FrontierEmbedding::apply(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  const bool final) const
    {
      pimpl->embedding_score(edge, features);
    }

    void FrontierEmbedding::apply_coarse(state_ptr_type& state,
					 const state_ptr_set_type& states,
					 const edge_type& edge,
					 feature_set_type& features,
					 const bool final) const
    { }
    
    void FrontierEmbedding::apply_predict(state_ptr_type& state,
					  const state_ptr_set_type& states,
					  const edge_type& edge,
					  feature_set_type& features,
					  const bool final) const
    { }
    
    void FrontierEmbedding::apply_scan(state_ptr_type& state,
				       const state_ptr_set_type& states,
				       const edge_type& edge,
				       const int dot,
				       feature_set_type& features,
				       const bool final) const
    { }
    
    void FrontierEmbedding::apply_complete(state_ptr_type& state,
					   const state_ptr_set_type& states,
					   const edge_type& edge,
					   feature_set_type& features,
					   const bool final) const
    {
      apply(state, states, edge, features, final);
    }
  };
};
