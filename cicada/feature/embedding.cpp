//
//  Copyright(C) 2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <memory>

#include "cicada/feature/embedding.hpp"
#include "cicada/word_embedding.hpp"
#include "cicada/parameter.hpp"
#include "cicada/symbol_vector.hpp"

#include "utils/hashmurmur3.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"
#include "utils/random_seed.hpp"

namespace cicada
{
  namespace feature
  {
    class EmbeddingImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

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

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_set_type::feature_type feature_type;
      
      typedef std::vector<feature_type, std::allocator<feature_type> > feature_name_set_type;
      
      typedef rule_type::word_type word_type;
      
    public:
      EmbeddingImpl(const path_type& path, const std::string& name)
	: embedding_(&word_embedding_type::create(path)), path_(path)
      {
	initialize(name);
      }
      
      EmbeddingImpl(const EmbeddingImpl& x)
	: embedding_(&word_embedding_type::create(x.path_)),
	  path_(x.path_),
	  skip_sgml_tag_(x.skip_sgml_tag_),
	  feature_names_(x.feature_names_)
	  
      { }
      
      EmbeddingImpl& operator=(const EmbeddingImpl& x)
      {
	
	embedding_     = &word_embedding_type::create(x.path_);
	path_          = x.path_;
	skip_sgml_tag_ = x.skip_sgml_tag_;
	feature_names_ = x.feature_names_;
	
	return *this;
      }
      
      void initialize(const std::string& name)
      {
	if (! embedding_)
	  throw std::runtime_error("no embedding?");
	
	feature_names_.clear();
	for (size_type i = 0; i != embedding_->dimension(); ++ i)
	  feature_names_.push_back(name + ':' + utils::lexical_cast<std::string>(i));
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
	  embedding_score(edge.rule->rhs.begin(), edge.rule->rhs.end(), features, extract_word(), skipper_sgml());
	else
	  embedding_score(edge.rule->rhs.begin(), edge.rule->rhs.end(), features, extract_word(), skipper_epsilon());
      }

      
      template <typename Iterator, typename Extract, typename Skipper>
      void embedding_score(Iterator first, Iterator last, 
			   feature_set_type& features,
			   Extract extract,
			   Skipper skipper) const
      {
	tensor_type layer = tensor_type::Zero(embedding_->dimension(), 1);
	
	for (/**/; first != last; ++ first)
	  if (first->is_terminal() && ! skipper(*first))
	    layer += embedding_->operator()(extract(*first));
	
	for (size_type i = 0; i != embedding_->dimension(); ++ i)
	  if (layer(i, 0) != parameter_type(0))
	    features[feature_names_[i]] = layer(i, 0);
	  else
	    features.erase(feature_names_[i]);
      }
      
    public:
      const word_embedding_type* embedding_;
      
      path_type path_;

      bool skip_sgml_tag_;
      
      // names...
      feature_name_set_type feature_names_;
    };
    
    Embedding::Embedding(const std::string& parameter)
      : pimpl(0)
    {
      typedef boost::filesystem::path path_type;
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "embedding")
	throw std::runtime_error("is this really embedding feature function? " + parameter);
      
      path_type   path;
      bool        skip_sgml_tag = false;
      std::string name;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for ngram: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (! boost::filesystem::exists(path))
	throw std::runtime_error("no embedding file? " + path.string());
      
      if (name.empty())
	name = "embedding";
      
      std::unique_ptr<impl_type> embedding_impl(new impl_type(path, name));
      
      embedding_impl->skip_sgml_tag_ = skip_sgml_tag;
      
      base_type::__state_size   = 0;
      base_type::__feature_name = name;
      
      pimpl = embedding_impl.release();
    }
    
    Embedding::~Embedding()
    {
      std::unique_ptr<impl_type> tmp(pimpl);
    }
    
    Embedding::Embedding(const Embedding& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}
    
    Embedding& Embedding::operator=(const Embedding& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    
    void Embedding::apply(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  const bool final) const
    {
      pimpl->embedding_score(edge, features);
    }

    void Embedding::apply_coarse(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 const bool final) const
    { }
    
    void Embedding::apply_predict(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  const bool final) const
    { }
    
    void Embedding::apply_scan(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       const int dot,
			       feature_set_type& features,
			       const bool final) const
    { }
    
    void Embedding::apply_complete(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   const bool final) const
    {
      apply(state, states, edge, features, final);
    }
  };
};
