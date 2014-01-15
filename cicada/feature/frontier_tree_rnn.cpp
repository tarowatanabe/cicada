//
//  Copyright(C) 2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <memory>

#include "cicada/feature/frontier_tree_rnn.hpp"
#include "cicada/bitree_rnn.hpp"
#include "cicada/parameter.hpp"
#include "cicada/symbol_vector.hpp"

#include "utils/array_power2.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"
#include "utils/small_vector.hpp"
#include "utils/random_seed.hpp"
#include "utils/space_separator.hpp"
#include "utils/array_power2.hpp"

#include <boost/tokenizer.hpp>

namespace cicada
{
  namespace feature
  {
    class FrontierTreeRNNImpl : public utils::hashmurmur3<size_t>
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef utils::hashmurmur3<size_t> hasher_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Vocab  vocab_type;
      
      typedef cicada::BiTreeRNN rnn_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;

      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;

      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
      typedef cicada::feature::FrontierTreeRNN::feature_name_set_type feature_name_set_type;

      typedef rnn_type::parameter_type parameter_type;
      typedef rnn_type::matrix_type    matrix_type;
      typedef rnn_type::tensor_type    tensor_type;

      typedef rnn_type::path_type path_type;

      typedef rule_type::word_type word_type;
      
      typedef std::vector<parameter_type, std::allocator<parameter_type> > buffer_type;
      typedef std::vector<word_type, std::allocator<word_type> > word_set_type;

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
      FrontierTreeRNNImpl(const path_type& path, const std::string& name)
	: rnn(path),
	  no_bos_eos(false), skip_sgml_tag(false),
	  attr_frontier_source("frontier-source"),
	  attr_frontier_target("frontier-target")
      {
	initialize(name);
      }
      
      FrontierTreeRNNImpl(const size_type& hidden,
			  const size_type& embedding,
			  const path_type& path_source,
			  const path_type& path_target,
			  const std::string& name)
	: rnn(hidden, embedding, path_source, path_target),
	  no_bos_eos(false), skip_sgml_tag(false),
	  attr_frontier_source("frontier-source"),
	  attr_frontier_target("frontier-target")
      {
	initialize(name);
      }
      
      FrontierTreeRNNImpl(const FrontierTreeRNNImpl& x)
	: rnn(x.rnn),
	  no_bos_eos(x.no_bos_eos),
	  skip_sgml_tag(x.skip_sgml_tag),
	  feature_names(x.feature_names),
	  attr_frontier_source("frontier-source"),
	  attr_frontier_target("frontier-target")
      {
      }

      FrontierTreeRNNImpl& operator=(const FrontierTreeRNNImpl& x)
      {
	rnn = x.rnn;
	
	no_bos_eos    = x.no_bos_eos;
	skip_sgml_tag = x.skip_sgml_tag;
	feature_names = x.feature_names;
	
	return *this;
      }
      
      void initialize(const std::string& name)
      {
	feature_names.clear();
	for (size_type i = 0; i != rnn.hidden_; ++ i)
	  feature_names.push_back(name + ':' + utils::lexical_cast<std::string>(i));
      }

      void initialize()
      {
	// pre-apply f(x) for Bi
	init = rnn.Bi_.array().unaryExpr(rnn_type::shtanh());
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
      
      void rnn_score(state_ptr_type& state,
		     const state_ptr_set_type& states,
		     const edge_type& edge,
		     feature_set_type& features,
		     const bool final) const
      {
	if (final && ! no_bos_eos) {
	  word_set_type& words = const_cast<word_set_type&>(words_tmp);
	  
	  words.clear();
	  words.push_back(vocab_type::BOS);
	  words.insert(words.end(), edge.rule->rhs.begin(), edge.rule->rhs.end());
	  words.push_back(vocab_type::EOS);
	  
	  if (skip_sgml_tag)
	    rnn_score(state, states, edge, words.begin(), words.end(),  features, extract_word(), skipper_sgml());
	  else
	    rnn_score(state, states, edge, words.begin(), words.end(), features, extract_word(), skipper_epsilon());
	} else {
	  if (skip_sgml_tag)
	    rnn_score(state, states, edge, edge.rule->rhs.begin(), edge.rule->rhs.end(), features, extract_word(), skipper_sgml());
	  else
	    rnn_score(state, states, edge, edge.rule->rhs.begin(), edge.rule->rhs.end(), features, extract_word(), skipper_epsilon());
	}
      }

      
      template <typename Iterator, typename Extract, typename Skipper>
      void rnn_score(state_ptr_type& state,
		     const state_ptr_set_type& states,
		     const edge_type& edge,
		     Iterator first, Iterator last, 
		     feature_set_type& features,
		     Extract extract,
		     Skipper skipper) const
      {
	const size_type offset1 = 0;
	const size_type offset2 = rnn.hidden_;
	const size_type offset_source = rnn.hidden_;
	const size_type offset_target = rnn.hidden_ + rnn.embedding_;
	
	const_cast<buffer_type&>(buffer_tmp).resize(rnn.hidden_);
	
	parameter_type* pointer_curr = const_cast<parameter_type*>(&(*buffer_tmp.begin()));
	parameter_type* pointer_next = reinterpret_cast<parameter_type*>(state);
	
	bool is_initial = true;
	int non_terminal_pos = 0;
	for (/**/; first != last; ++ first)
	  if (first->is_non_terminal()) {
	    matrix_type buffer_curr(pointer_curr, rnn.hidden_, 1);
	    matrix_type buffer_next(pointer_next, rnn.hidden_, 1);
	    
	    const int __non_terminal_index = first->non_terminal_index();
	    const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0,
								non_terminal_pos,
								__non_terminal_index - 1);
	    ++ non_terminal_pos;
	    
	    matrix_type buffer_ante(const_cast<parameter_type*>(reinterpret_cast<const parameter_type*>(states[antecedent_index])),
				    rnn.hidden_, 1);
	    
	    if (is_initial)
	      buffer_next = (rnn.Bn_
			     + rnn.Wn_.block(0, offset1, rnn.hidden_, rnn.hidden_) * init
			     + rnn.Wn_.block(0, offset2, rnn.hidden_, rnn.hidden_) * buffer_ante
			     ).array().unaryExpr(rnn_type::shtanh());
	    else 
	      buffer_next = (rnn.Bn_
			     + rnn.Wn_.block(0, offset1, rnn.hidden_, rnn.hidden_) * buffer_curr
			     + rnn.Wn_.block(0, offset2, rnn.hidden_, rnn.hidden_) * buffer_ante
			     ).array().unaryExpr(rnn_type::shtanh());
	    
	    std::swap(pointer_curr, pointer_next);
	    is_initial = false;
	  } 
	
	// iterate over phrase-pairs...
	attribute_set_type::const_iterator siter = edge.attributes.find(attr_frontier_source);
	attribute_set_type::const_iterator titer = edge.attributes.find(attr_frontier_target);
	
	if (siter != edge.attributes.end() && titer != edge.attributes.end()) {
	  const std::string& frontier_source = boost::apply_visitor(__attribute_string(), siter->second);
	  const std::string& frontier_target = boost::apply_visitor(__attribute_string(), titer->second);
	  
	  const phrase_type& phrase_source = cache_phrase(frontier_source, cache_source, skipper);
	  const phrase_type& phrase_target = cache_phrase(frontier_target, cache_target, skipper);
	  
	  if (! phrase_source.empty() || ! phrase_target.empty()) {
	    
	    if (phrase_source.empty()) {
	      // only target...
	      phrase_type::const_iterator titer_begin = phrase_target.begin();
	      phrase_type::const_iterator titer_end   = phrase_target.end();
	      for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
		matrix_type buffer_curr(pointer_curr, rnn.hidden_, 1);
		matrix_type buffer_next(pointer_next, rnn.hidden_, 1);
		
		if (is_initial)
		  buffer_next = (rnn.Bt_
				 + rnn.Wt_.block(0, offset1, rnn.hidden_, rnn.hidden_)    * init
				 + rnn.Wt_.block(0, offset2, rnn.hidden_, rnn.embedding_) * rnn.target_->operator()(*titer)
				 ).array().unaryExpr(rnn_type::shtanh());
		else
		  buffer_next = (rnn.Bt_
				 + rnn.Wt_.block(0, offset1, rnn.hidden_, rnn.hidden_)    * buffer_curr
				 + rnn.Wt_.block(0, offset2, rnn.hidden_, rnn.embedding_) * rnn.target_->operator()(*titer)
				 ).array().unaryExpr(rnn_type::shtanh());
		
		std::swap(pointer_curr, pointer_next);
		is_initial = false;
	      }
	      
	    } else if (phrase_target.empty()) {
	      // only source...
	      phrase_type::const_iterator siter_begin = phrase_source.begin();
	      phrase_type::const_iterator siter_end   = phrase_source.end();
	      for (phrase_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
		matrix_type buffer_curr(pointer_curr, rnn.hidden_, 1);
		matrix_type buffer_next(pointer_next, rnn.hidden_, 1);

		if (is_initial)
		  buffer_next = (rnn.Bs_
				 + rnn.Ws_.block(0, offset1, rnn.hidden_, rnn.hidden_)    * init
				 + rnn.Ws_.block(0, offset2, rnn.hidden_, rnn.embedding_) * rnn.source_->operator()(*siter)
				 ).array().unaryExpr(rnn_type::shtanh());
		else
		  buffer_next = (rnn.Bs_
				 + rnn.Ws_.block(0, offset1, rnn.hidden_, rnn.hidden_)    * buffer_curr
				 + rnn.Ws_.block(0, offset2, rnn.hidden_, rnn.embedding_) * rnn.source_->operator()(*siter)
				 ).array().unaryExpr(rnn_type::shtanh());
		
		std::swap(pointer_curr, pointer_next);
		is_initial = false;
	      }
	    } else {
	      // pairs...
	      phrase_type::const_iterator siter_begin = phrase_source.begin();
	      phrase_type::const_iterator siter_end   = phrase_source.end();
	      phrase_type::const_iterator titer_begin = phrase_target.begin();
	      phrase_type::const_iterator titer_end   = phrase_target.end();
	      
	      for (phrase_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
		for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
		  matrix_type buffer_curr(pointer_curr, rnn.hidden_, 1);
		  matrix_type buffer_next(pointer_next, rnn.hidden_, 1);
		  
		  if (is_initial)
		    buffer_next = (rnn.Bp_
				   + rnn.Wp_.block(0, offset1, rnn.hidden_, rnn.hidden_)          * init
				   + rnn.Wp_.block(0, offset_source, rnn.hidden_, rnn.embedding_) * rnn.source_->operator()(*siter)
				   + rnn.Wp_.block(0, offset_target, rnn.hidden_, rnn.embedding_) * rnn.target_->operator()(*titer)
				   ).array().unaryExpr(rnn_type::shtanh());
		  else
		    buffer_next = (rnn.Bp_
				   + rnn.Wp_.block(0, offset1, rnn.hidden_, rnn.hidden_)          * buffer_curr
				   + rnn.Wp_.block(0, offset_source, rnn.hidden_, rnn.embedding_) * rnn.source_->operator()(*siter)
				   + rnn.Wp_.block(0, offset_target, rnn.hidden_, rnn.embedding_) * rnn.target_->operator()(*titer)
				   ).array().unaryExpr(rnn_type::shtanh());
		  
		  std::swap(pointer_curr, pointer_next);
		  is_initial = false;
		}
	    }
	  }
	}
	
	// copy into state buffer when necessary..
	if (is_initial) {
	  // nothing is propagated... this may not happen, but we will simply copy "init"
	  matrix_type buffer_next(reinterpret_cast<parameter_type*>(state), rnn.hidden_, 1);
	  
	  buffer_next = init;
	} else if (pointer_curr != reinterpret_cast<parameter_type*>(state)) {
	  matrix_type buffer_curr(pointer_curr, rnn.hidden_, 1);
	  matrix_type buffer_next(reinterpret_cast<parameter_type*>(state), rnn.hidden_, 1);
	  
	  buffer_next = buffer_curr;
	}
	
	// add features...
	features.reserve(features.size() + rnn.hidden_);
	
	matrix_type buffer(reinterpret_cast<parameter_type*>(state), rnn.hidden_, 1);
	for (size_type i = 0; i != rnn.hidden_; ++ i) {
	  if (buffer(i, 0) != parameter_type(0))
	    features[feature_names[i]] = buffer(i, 0);
	  else
	    features.erase(feature_names[i]);
	}
      }

    private:
      template <typename Skipper>
      const phrase_type& cache_phrase(const std::string& frontier,
				      const cache_phrase_set_type& caches,
				      Skipper skipper) const
      {
	if (frontier.empty()) return phrase_tmp;
	
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
      rnn_type    rnn;
      tensor_type init;

      buffer_type   buffer_tmp;
      word_set_type words_tmp;
      phrase_type   phrase_tmp;

      cache_phrase_set_type cache_source;
      cache_phrase_set_type cache_target;
      
      bool no_bos_eos;
      bool skip_sgml_tag;
      
      // names...
      feature_name_set_type feature_names;

      attribute_type attr_frontier_source;
      attribute_type attr_frontier_target;
    };
    
    FrontierTreeRNN::FrontierTreeRNN(const std::string& parameter)
      : pimpl(0)
    {
      typedef boost::filesystem::path path_type;
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "frontier-tree-rnn")
	throw std::runtime_error("is this really frontier tree-rnn feature function? " + parameter);
      
      path_type   path;
      path_type   path_source;
      path_type   path_target;
      size_type   hidden = 0;
      size_type   embedding = 0;
      bool        skip_sgml_tag = false;
      bool        no_bos_eos = false;
      
      std::string name;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else if (utils::ipiece(piter->first) == "embedding-source-file")
	  path_source = piter->second;
	else if (utils::ipiece(piter->first) == "embedding-target-file")
	  path_target = piter->second;
	else if (utils::ipiece(piter->first) == "dimension-hidden")
	  hidden = utils::lexical_cast<size_type>(piter->second);
	else if (utils::ipiece(piter->first) == "dimension-embedding")
	  embedding = utils::lexical_cast<size_type>(piter->second);
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "no-bos-eos")
	  no_bos_eos = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for ngram: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (! path.empty() && ! path_source.empty() && ! path_target.empty())
	throw std::runtime_error("either one of model file or source/target embedding file");
      
      if (! path.empty()) {
	if (! boost::filesystem::exists(path))
	  throw std::runtime_error("no model file? " + path.string());
	
      } else if (! path_source.empty() && ! path_target.empty()) {
	if (! boost::filesystem::exists(path_source))
	  throw std::runtime_error("no source embedding file? " + path_source.string());
	if (! boost::filesystem::exists(path_target))
	  throw std::runtime_error("no target embedding file? " + path_target.string());
	
	if (hidden == 0 || embedding == 0)
	  throw std::runtime_error("invalid dimension");
	
      } else
	throw std::runtime_error("no model file nor source/target embedding file?");
      
      if (name.empty())
	name = "frontier-tree-rnn";

      const bool model_mode = (! path.empty());
      
      std::auto_ptr<impl_type> rnn_impl(model_mode
					? new impl_type(path, name)
					: new impl_type(hidden, embedding, path_source, path_target, name));
      
      rnn_impl->no_bos_eos    = no_bos_eos;
      rnn_impl->skip_sgml_tag = skip_sgml_tag;
      
      if (! model_mode) {
	boost::mt19937 generator;
	generator.seed(utils::random_seed());
	
	rnn_impl->rnn.random(generator);
      }
      
      base_type::__state_size = sizeof(tree_rnn_type::parameter_type) * rnn_impl->rnn.hidden_;
      base_type::__feature_name = name;
      
      pimpl = rnn_impl.release();
    }
    
    FrontierTreeRNN::~FrontierTreeRNN()
    {
      std::auto_ptr<impl_type> tmp(pimpl);
    }
    
    FrontierTreeRNN::FrontierTreeRNN(const FrontierTreeRNN& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    FrontierTreeRNN::tree_rnn_type& FrontierTreeRNN::model() const
    {
      return const_cast<tree_rnn_type&>(pimpl->rnn);
    }

    const FrontierTreeRNN::feature_name_set_type& FrontierTreeRNN::features() const
    {
      return pimpl->feature_names;
    }
    
    bool FrontierTreeRNN::no_bos_eos() const
    {
      return pimpl->no_bos_eos;
    }

    bool FrontierTreeRNN::skip_sgml_tag() const
    {
      return pimpl->skip_sgml_tag;
    }

    FrontierTreeRNN& FrontierTreeRNN::operator=(const FrontierTreeRNN& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void FrontierTreeRNN::initialize()
    {
      pimpl->initialize();
    }
    
    void FrontierTreeRNN::apply(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				const bool final) const
    {
      pimpl->rnn_score(state, states, edge, features, final);
    }

    void FrontierTreeRNN::apply_coarse(state_ptr_type& state,
				       const state_ptr_set_type& states,
				       const edge_type& edge,
				       feature_set_type& features,
				       const bool final) const
    { }
    
    void FrontierTreeRNN::apply_predict(state_ptr_type& state,
					const state_ptr_set_type& states,
					const edge_type& edge,
					feature_set_type& features,
					const bool final) const
    { }
    
    void FrontierTreeRNN::apply_scan(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     const int dot,
				     feature_set_type& features,
				     const bool final) const
    { }
    
    void FrontierTreeRNN::apply_complete(state_ptr_type& state,
					 const state_ptr_set_type& states,
					 const edge_type& edge,
					 feature_set_type& features,
					 const bool final) const
    {
      apply(state, states, edge, features, final);
    }
  };
};
