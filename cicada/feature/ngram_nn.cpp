//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <memory>

#include "cicada/feature/ngram_nn.hpp"
#include "cicada/ngram_nn.hpp"
#include "cicada/ngram_cache.hpp"
#include "cicada/parameter.hpp"
#include "cicada/symbol_vector.hpp"

#include "utils/array_power2.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"
#include "utils/small_vector.hpp"

namespace cicada
{
  namespace feature
  {
    class NGramNNImpl : public utils::hashmurmur3<size_t>
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Vocab  vocab_type;
      
      typedef cicada::NGramNN      ngram_type;
      
      typedef cicada::NGramCache<symbol_type::id_type, double> ngram_cache_type;
      
      typedef std::vector<symbol_type::id_type, std::allocator<symbol_type::id_type> > buffer_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_set_type::feature_type feature_type;
      
      typedef rule_type::symbol_set_type phrase_type;
      
      typedef utils::hashmurmur3<size_t> hasher_type;


      struct CacheContext
      {
	typedef utils::small_vector<symbol_type::id_type, std::allocator<symbol_type::id_type> > phrase_type;
	
	phrase_type context;
	phrase_type ngram;
	double      score;
	
	CacheContext() : context(), ngram(), score() {}
      };
      
      typedef CacheContext cache_context_type;
      typedef utils::array_power2<cache_context_type,  1024 * 128, std::allocator<cache_context_type> >  cache_context_set_type;
      
      typedef boost::filesystem::path path_type;

      struct result_type
      {
	result_type() : prob_(0), bound_(0), oov_(0) {}
	
	double prob_;
	double bound_;
	int    oov_;
      };
      
    public:
      // we will use -2, since -1 may be used for OOV
      static const symbol_type::id_type id_empty;
      static const symbol_type::id_type id_star;

    public:
      NGramNNImpl(const path_type& __path, const int __order, const bool populate)
	: ngram(&ngram_type::create(__path)),
	  order(__order), no_bos_eos(false), skip_sgml_tag(false), split_estimate(false)
      {
	order = utils::bithack::min(order, ngram->order());

	if (populate)
	  ngram->populate();
	
	initialize_cache();
	
	id_oov = ngram->vocab()[vocab_type::UNK];
	id_bos = ngram->vocab()[vocab_type::BOS];
	id_eos = ngram->vocab()[vocab_type::EOS];
      }

      NGramNNImpl(const NGramNNImpl& x)
	: ngram(&ngram_type::create(x.ngram->path())),
	  order(x.order),
	  no_bos_eos(x.no_bos_eos),
	  skip_sgml_tag(x.skip_sgml_tag),
	  split_estimate(x.split_estimate),
	  feature_name(x.feature_name),
	  feature_name_oov(x.feature_name_oov),
	  feature_name_estimate(x.feature_name_estimate),
	  id_oov(x.id_oov),
	  id_bos(x.id_bos),
	  id_eos(x.id_eos)
	  
      {
	initialize_cache();
      }

      NGramNNImpl& operator=(const NGramNNImpl& x)
      {
	ngram = &ngram_type::create(x.ngram->path());
	order = x.order;
	
	no_bos_eos     = x.no_bos_eos;
	skip_sgml_tag  = x.skip_sgml_tag;
	split_estimate = x.split_estimate;
	
	feature_name          = x.feature_name;
	feature_name_oov      = x.feature_name_oov;
	feature_name_estimate = x.feature_name_estimate;
	
	id_oov           = x.id_oov;
	id_bos           = x.id_bos;
	id_eos           = x.id_eos;
	
	initialize_cache();
	
	return *this;
      }
      
      void initialize_cache()
      {
	cache_logprob.clear();
	cache_estimate = ngram_cache_type(ngram->order());
      }
      
      template <typename Iterator>
      inline
      size_t hash_phrase(Iterator first, Iterator last, size_t seed=0) const
      {
	return hasher_type::operator()(first, last, seed);
      }
      
      template <typename Iterator, typename __Phrase>
      inline
      bool equal_phrase(Iterator first, Iterator last, const __Phrase& x) const
      {
	return static_cast<int>(x.size()) == std::distance(first, last) && std::equal(first, last, x.begin());
      }
      
      template <typename Iterator>
      double ngram_score(Iterator first, Iterator iter, Iterator last) const
      {
	if (iter == last) return 0.0;
	
	first = std::max(first, iter - order + 1);
	
	const size_t cache_pos = hash_phrase(first, iter, hash_phrase(iter, last)) & (cache_logprob.size() - 1);
	cache_context_type& cache = const_cast<cache_context_type&>(cache_logprob[cache_pos]);
	
	if (! equal_phrase(first, iter, cache.context) || ! equal_phrase(iter, last, cache.ngram)) {
	  cache.context.assign(first, iter);
	  cache.ngram.assign(iter, last);
	  
	  buffer_type& buffer = const_cast<buffer_type&>(buffer_score_impl);
	  buffer.clear();
	  buffer.insert(buffer.end(), first, iter);
	  
	  cache.score = 0.0;
	  for (/**/; iter != last; ++ iter) {
	    buffer.push_back(*iter);
	    
	    cache.score += ngram->operator()(std::max(buffer.begin(), buffer.end() - order), buffer.end());
	  }
	}
	
	return cache.score;
      }
      
      template <typename Iterator>
      double ngram_estimate(Iterator first, Iterator last) const
      {
	if (first == last) return 0.0;
	
	const size_type cache_pos = cache_estimate(first, last);
	
	if (! cache_estimate.equal_to(cache_pos, first, last)) {
	  ngram_cache_type& cache = const_cast<ngram_cache_type&>(cache_estimate);

	  cache.assign(cache_pos, first, last);
	  
	  buffer_type& buffer = const_cast<buffer_type&>(buffer_score_impl);
	  buffer.clear();
	  
	  // skip initial bos...
	  if (*first == id_bos) {
	    buffer.push_back(*first);
	    ++ first;
	  }
	  
	  double score = 0.0;
	  for (/**/; first != last; ++ first) {
	    buffer.push_back(*first);
	    
	    score += ngram->operator()(buffer.begin(), buffer.end());
	  }
	  
	  cache[cache_pos] = score;
	}
	
	return cache_estimate[cache_pos];
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
      
      void ngram_score(state_ptr_type& state,
		       const state_ptr_set_type& states,
		       const edge_type& edge,
		       result_type& result) const
      {
	if (skip_sgml_tag)
	  ngram_score(state, states, edge, result, extract_word(), skipper_sgml());
	else
	  ngram_score(state, states, edge, result, extract_word(), skipper_epsilon());
      }

      void ngram_estimate(const edge_type& edge, result_type& result) const
      {
	if (skip_sgml_tag)
	  ngram_estimate(edge, result, extract_word(), skipper_sgml());
	else
	  ngram_estimate(edge, result, extract_word(), skipper_epsilon());
      }
      
      template <typename Extract, typename Skipper>
      void ngram_score(state_ptr_type& state,
		       const state_ptr_set_type& states,
		       const edge_type& edge,
		       result_type& result,
		       Extract extract,
		       Skipper skipper) const
      {
	const int context_size = order - 1;
	const rule_type& rule = *(edge.rule);
	const phrase_type& target = rule.rhs;
	
	phrase_type::const_iterator titer_begin = target.begin();
	phrase_type::const_iterator titer_end   = target.end();
	
	// we will reserve enough space so that buffer's memory will not be re-allocated.
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.reserve((titer_end - titer_begin) + (order * 2) * states.size());
	buffer.clear();
	
	if (states.empty()) {
	  // we will copy to buffer...
	  for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	    if (! skipper(*titer)) {
	      buffer.push_back(ngram->vocab()[extract(*titer)]);
	      result.oov_ += (buffer.back() == id_oov);
	    }
	  
	  symbol_type::id_type* context = reinterpret_cast<symbol_type::id_type*>(state);
	  symbol_type::id_type* context_end = context + order * 2;
	  
	  if (static_cast<int>(buffer.size()) <= context_size) {
	    std::fill(std::copy(buffer.begin(), buffer.end(), context), context_end, id_empty);
	    
	    result.bound_ += ngram_estimate(buffer.begin(), buffer.end());
	  } else {
	    buffer_type::const_iterator biter_begin = buffer.begin();
	    buffer_type::const_iterator biter_end   = buffer.end();
	    
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix(biter_begin, biter_begin + context_size);
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> suffix(biter_end - context_size, biter_end);
	    
	    std::copy(prefix.first, prefix.second, context);
	    context[prefix.second - prefix.first] = id_star;
	    std::fill(std::copy(suffix.first, suffix.second, context + (prefix.second - prefix.first) + 1), context_end, id_empty);
	    
	    result.bound_ += ngram_estimate(prefix.first, prefix.second);
	    result.prob_  += ngram_score(prefix.first, prefix.second, biter_end);
	  }
	} else {
	  int star_first = -1;
	  int star_last  = -1;
	  
	  buffer_type::iterator biter_first = buffer.begin();
	  buffer_type::iterator biter       = buffer.begin();
	  
	  int non_terminal_pos = 0;
	  for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	    if (titer->is_non_terminal()) {
	      const int __non_terminal_index = titer->non_terminal_index();
	      const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	      ++ non_terminal_pos;
	      
	      const symbol_type::id_type* context = reinterpret_cast<const symbol_type::id_type*>(states[antecedent_index]);
	      const symbol_type::id_type* context_end  = std::find(context, context + order * 2, id_empty);
	      const symbol_type::id_type* context_star = std::find(context, context_end, id_star);

	      // adjust antecedent context
	      result.bound_ -= ngram_estimate(context, context_star);
	      
	      buffer.insert(buffer.end(), context, context_star);
	      if (biter - biter_first >= context_size || star_first >= 0)
		result.prob_ += ngram_score(biter_first, biter, buffer.end());
	      else
		result.prob_ += ngram_score(biter_first, std::min(biter_first + context_size, buffer.end()), buffer.end());
	      biter = buffer.end();
	      
	      if (context_star != context_end) {
		star_first = utils::bithack::branch(star_first < 0, static_cast<int>(buffer.size()), star_first);
		star_last  = buffer.size();
		
		biter_first = buffer.end() + 1;
		buffer.insert(buffer.end(), context_star, context_end);
		biter = buffer.end();
	      }
	      
	    } else if (! skipper(*titer)) {
	      buffer.push_back(ngram->vocab()[extract(*titer)]);
	      result.oov_ += (buffer.back() == id_oov);
	    }
	  }
	  
	  if (biter != buffer.end()) {
	    if (biter - biter_first >= context_size || star_first >= 0)
	      result.prob_ += ngram_score(biter_first, biter, buffer.end());
	    else
	      result.prob_ += ngram_score(biter_first, std::min(biter_first + context_size, buffer.end()), buffer.end());
	  }
	  
	  // construct state vector..
	  symbol_type::id_type* context = reinterpret_cast<symbol_type::id_type*>(state);
	  symbol_type::id_type* context_end = context + order * 2;
	  
	  if (star_first >= 0) {
	    if (buffer[star_first] != id_star)
	      throw std::runtime_error("no star at star-first?");
	    if (buffer[star_last] != id_star)
	      throw std::runtime_error("no star at star-last?");

	    const int prefix_size = utils::bithack::min(star_first, context_size);
	    const int suffix_size = utils::bithack::min(int(buffer.size() - (star_last + 1)), context_size);
	    
	    buffer_type::const_iterator biter_begin = buffer.begin();
	    buffer_type::const_iterator biter_end   = buffer.end();
	    
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix(biter_begin, biter_begin + prefix_size);
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> suffix(biter_end - suffix_size, biter_end);
	    
	    std::copy(prefix.first, prefix.second, context);
	    context[prefix.second - prefix.first] = id_star;
	    std::fill(std::copy(suffix.first, suffix.second, context + (prefix.second - prefix.first) + 1), context_end, id_empty);
	    
	    result.bound_ += ngram_estimate(prefix.first, prefix.second);
	    result.prob_  += ngram_score(prefix.first, prefix.second, biter_begin + prefix_size);
	  } else {
	    if (static_cast<int>(buffer.size()) <= context_size) {
	      std::fill(std::copy(buffer.begin(), buffer.end(), context), context_end, id_empty);
	      
	      result.bound_ += ngram_estimate(buffer.begin(), buffer.end());
	    } else {
	      buffer_type::const_iterator biter_begin = buffer.begin();
	      buffer_type::const_iterator biter_end   = buffer.end();
	      
	      std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix(biter_begin, biter_begin + context_size);
	      std::pair<buffer_type::const_iterator, buffer_type::const_iterator> suffix(biter_end - context_size, biter_end);
	      
	      std::copy(prefix.first, prefix.second, context);
	      context[prefix.second - prefix.first] = id_star;
	      std::fill(std::copy(suffix.first, suffix.second, context + (prefix.second - prefix.first) + 1), context_end, id_empty);
	      
	      result.bound_ += ngram_estimate(prefix.first, prefix.second);
	      result.prob_  += ngram_score(prefix.first, prefix.second, biter_begin + context_size);
	    }
	  }
	}
      }
      
      template <typename Extract, typename Skipper>
      void ngram_estimate(const edge_type& edge, result_type& result, Extract extract, Skipper skipper) const
      {
	const int context_size = order - 1;
	const rule_type& rule = *(edge.rule);
	const phrase_type& target = rule.rhs;
	
	phrase_type::const_iterator titer_begin = target.begin();
	phrase_type::const_iterator titer_end   = target.end();
	
	// we will reserve enough space so that buffer's memory will not be re-allocated.
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	buffer.reserve(titer_end - titer_begin);
	
	for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	  if (titer->is_non_terminal()) {
	    if (! buffer.empty()) {
	      buffer_type::iterator biter_begin = buffer.begin();
	      buffer_type::iterator biter_end   = buffer.end();
	      buffer_type::iterator biter       = std::min(biter_begin + context_size, biter_end);
	      
	      result.bound_ += ngram_estimate(biter_begin, biter);
	      result.prob_  += ngram_score(biter_begin, biter, biter_end);
	    }
	    
	    buffer.clear();
	  } else if (! skipper(*titer)) {
	    buffer.push_back(ngram->vocab()[extract(*titer)]);
	    result.oov_ += (buffer.back() == id_oov);
	  }
	}
	
	if (! buffer.empty()) {
	  buffer_type::iterator biter_begin = buffer.begin();
	  buffer_type::iterator biter_end   = buffer.end();
	  buffer_type::iterator biter       = std::min(biter_begin + context_size, biter_end);
	  
	  result.bound_ += ngram_estimate(biter_begin, biter);
	  result.prob_  += ngram_score(biter_begin, biter, biter_end);
	}
      }
      
      void ngram_final_score(const state_ptr_type& state, result_type& result) const
      {
	const symbol_type::id_type* context      = reinterpret_cast<const symbol_type::id_type*>(state);
	const symbol_type::id_type* context_end  = std::find(context, context + order * 2, id_empty);
	const symbol_type::id_type* context_star = std::find(context, context_end, id_star);
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	
	if (no_bos_eos) {
	  buffer.insert(buffer.end(), context, context_star);
	  
	  if (! buffer.empty())
	    result.prob_ += ngram_score(buffer.begin(), buffer.begin() + (buffer.front() == id_bos), buffer.end());
	  
	  result.bound_ -= ngram_estimate(buffer.begin(), buffer.end());
	} else {
	  buffer.push_back(id_bos);
	  buffer.insert(buffer.end(), context, context_star);

	  result.prob_  += ngram_score(buffer.begin(), buffer.begin() + 1, buffer.end());
	  result.bound_ -= ngram_estimate(buffer.begin() + 1, buffer.end());
	  
	  if (context_star != context_end) {
	    buffer.clear();
	    buffer.insert(buffer.end(), context_star + 1, context_end);
	  }
	  
	  buffer.push_back(id_eos);
	  
	  result.prob_ += ngram_score(buffer.begin(), buffer.end() - 1, buffer.end());
	}
      }
            
      void ngram_predict_score(const state_ptr_type& state, result_type& result)
      {
	symbol_type::id_type* context = reinterpret_cast<symbol_type::id_type*>(state);
	
	std::fill(context, context + order * 2, id_empty);
	
	// for no-bos-eos, we need to keep track of whether P(bos) was scored or not...
	if (! no_bos_eos)
	  context[0] = id_bos;
      }

      void ngram_scan_score(state_ptr_type& state,
			    const edge_type& edge,
			    const int dot,
			    result_type& result)
      {
	if (skip_sgml_tag)
	  ngram_scan_score(state, edge, dot, result, extract_word(), skipper_sgml());
	else
	  ngram_scan_score(state, edge, dot, result, extract_word(), skipper_epsilon());
      }
      
      template <typename Extract, typename Skipper>
      void ngram_scan_score(state_ptr_type& state,
			    const edge_type& edge,
			    const int dot,
			    result_type& result,
			    Extract extract,
			    Skipper skipper)
      {
	const int context_size = order - 1;

	const rule_type& rule = *(edge.rule);
	const phrase_type& phrase = rule.rhs;

	symbol_type::id_type* context = reinterpret_cast<symbol_type::id_type*>(state);
	symbol_type::id_type* context_end  = std::find(context, context + order, id_empty);
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	buffer.reserve(order + phrase.size());
	buffer.insert(buffer.end(), context, context_end);
	
	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin() + dot; piter != piter_end && ! piter->is_non_terminal(); ++ piter)
	  if (! skipper(*piter)) {
	    buffer.push_back(ngram->vocab()[extract(*piter)]);
	    result.oov_ += (buffer.back() == id_oov);
	  }
	
	std::pair<buffer_type::iterator, buffer_type::iterator> suffix(std::max(buffer.begin(), buffer.end() - context_size), buffer.end());
	std::fill(std::copy(suffix.first, suffix.second, context), context + order * 2, id_empty);
	
	if (no_bos_eos && context == context_end && ! buffer.empty() && buffer.front() == id_bos)
	  result.prob_ += ngram_score(buffer.begin(), buffer.begin() + 1, buffer.end());
	else
	  result.prob_ += ngram_score(buffer.begin(), buffer.begin() + (context_end - context), buffer.end());
      }
      
      void ngram_complete_score(state_ptr_type& state, result_type& result)
      {
	if (no_bos_eos)
	  return;
	
	symbol_type::id_type* context = reinterpret_cast<symbol_type::id_type*>(state);
	symbol_type::id_type* context_end  = std::find(context, context + order, vocab_type::EMPTY);
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	buffer.insert(buffer.end(), context, context_end);
	buffer.push_back(id_eos);
	
	result.prob_ += ngram_score(buffer.begin(), buffer.end() - 1, buffer.end());
      }
      
      // caching...
      cache_context_set_type cache_logprob;
      ngram_cache_type       cache_estimate;
      
      // actual buffers...
      buffer_type    buffer_impl;
      buffer_type    buffer_score_impl;
      
      // ngrams
      ngram_type*     ngram;
      int             order;
      
      bool no_bos_eos;
      bool skip_sgml_tag;
      bool split_estimate;
      
      // names...
      feature_type feature_name;
      feature_type feature_name_oov;
      feature_type feature_name_estimate;
      
      symbol_type::id_type id_oov;
      symbol_type::id_type id_bos;
      symbol_type::id_type id_eos;
    };

    const NGramNNImpl::symbol_type::id_type NGramNNImpl::id_empty = NGramNNImpl::symbol_type::id_type(-2);
    const NGramNNImpl::symbol_type::id_type NGramNNImpl::id_star  = NGramNNImpl::symbol_type::id_type(-3);
    
    NGramNN::NGramNN(const std::string& parameter)
      : pimpl(0), pimpl_coarse(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (utils::ipiece(param.name()) != "ngram-nn")
	throw std::runtime_error("is this really ngram feature function? " + parameter);

      path_type   path;
      bool        populate = false;
      int         order = 3;
      bool        skip_sgml_tag = false;
      bool        split_estimate = false;
      bool        no_bos_eos = false;
      
      path_type   coarse_path;
      bool        coarse_populate;
      int         coarse_order = 0;
      
      std::string name;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else if (utils::ipiece(piter->first) == "populate")
	  populate = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "no-bos-eos")
	  no_bos_eos = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "split-estimate")
	  split_estimate = utils::lexical_cast<bool>(piter->second);	
	else if (utils::ipiece(piter->first) == "coarse-file")
	  coarse_path = piter->second;
	else if (utils::ipiece(piter->first) == "coarse-populate")
	  coarse_populate = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "coarse-order")
	  coarse_order = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for ngram: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (path.empty())
	throw std::runtime_error("no ngram file? " + path.string());
      
      if (order <= 0)
	throw std::runtime_error("invalid ngram order: " + utils::lexical_cast<std::string>(order));

      if (coarse_order > order)
	throw std::runtime_error("invalid coarse order: coarse-order <= order");
      if (! coarse_path.empty() && ! boost::filesystem::exists(coarse_path))
	throw std::runtime_error("no coarse ngram language model? " + coarse_path.string());
      
      std::unique_ptr<impl_type> ngram_impl(new impl_type(path, order, populate));
      
      ngram_impl->no_bos_eos     = no_bos_eos;
      ngram_impl->skip_sgml_tag  = skip_sgml_tag;
      ngram_impl->split_estimate = split_estimate;
      
      // two contexts (order - 1) for each edge, with two separator..
      base_type::__state_size = sizeof(symbol_type::id_type) * ngram_impl->order * 2;
      base_type::__feature_name = (name.empty() ? std::string("ngram-nn") : name);
      
      ngram_impl->feature_name     = base_type::__feature_name;
      ngram_impl->feature_name_oov = static_cast<const std::string&>(base_type::__feature_name) + ":oov-penalty";
      if (split_estimate)
	ngram_impl->feature_name_estimate = static_cast<const std::string&>(base_type::__feature_name) + ":estimate";
      
      pimpl = ngram_impl.release();
      
      // ...
      if (coarse_order > 0 || ! coarse_path.empty()) {
	
	if (coarse_order <= 0)
	  throw std::runtime_error("coarse order must be non-zero!");
	
	if (! coarse_path.empty()) {
	  std::unique_ptr<impl_type> ngram_impl(new impl_type(coarse_path, coarse_order, coarse_populate));

	  ngram_impl->no_bos_eos = no_bos_eos;
	  ngram_impl->skip_sgml_tag = skip_sgml_tag;
	  
	  pimpl_coarse = ngram_impl.release();
	} else {
	  std::unique_ptr<impl_type> ngram_impl(new impl_type(*pimpl));
	  ngram_impl->order = coarse_order;
	  
	  pimpl_coarse = ngram_impl.release();
	}
      }
    }
    
    NGramNN::~NGramNN()
    {
      std::unique_ptr<impl_type> tmp(pimpl);
      if (pimpl_coarse)
	std::unique_ptr<impl_type> tmp_coarse(pimpl_coarse);
    }
    
    NGramNN::NGramNN(const NGramNN& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl)),
	pimpl_coarse(x.pimpl_coarse ? new impl_type(*x.pimpl_coarse) : 0)
    {}
    
    NGramNN& NGramNN::operator=(const NGramNN& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      *pimpl = *x.pimpl;
      if (x.pimpl_coarse) {
	if (pimpl_coarse)
	  *pimpl_coarse = *x.pimpl_coarse;
	else
	  pimpl_coarse = new impl_type(*x.pimpl_coarse);
      } else {
	if (pimpl_coarse)
	  delete pimpl_coarse;
	pimpl_coarse = 0;
      }
      
      return *this;
    }
    
    
    void NGramNN::apply(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			const bool final) const
    {
      impl_type::result_type result;
      
      pimpl->ngram_score(state, states, edge, result);
      if (final)
	pimpl->ngram_final_score(state, result);
      
      if (pimpl->split_estimate) {
	if (result.prob_ != 0.0)
	  features[pimpl->feature_name] = result.prob_;
	else
	  features.erase(pimpl->feature_name);
	
	if (result.bound_ != 0.0)
	  features[pimpl->feature_name_estimate] = result.bound_;
	else
	  features.erase(pimpl->feature_name_estimate);
      } else {
	const double score = result.prob_ + result.bound_;
	
	if (score != 0.0)
	  features[pimpl->feature_name] = score;
	else
	  features.erase(pimpl->feature_name);
      }

      if (result.oov_)
	features[pimpl->feature_name_oov] = - result.oov_;
      else
	features.erase(pimpl->feature_name_oov);
    }

    void NGramNN::apply_coarse(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       const bool final) const
    {
      impl_type::result_type result;

      if (pimpl_coarse) {
	pimpl_coarse->ngram_score(state, states, edge, result);
	if (final)
	  pimpl_coarse->ngram_final_score(state, result);
      } else
	pimpl->ngram_estimate(edge, result);
      
      if (pimpl->split_estimate) {
	if (result.prob_ != 0.0)
	  features[pimpl->feature_name] = result.prob_;
	else
	  features.erase(pimpl->feature_name);
	
	if (result.bound_ != 0.0)
	  features[pimpl->feature_name_estimate] = result.bound_;
	else
	  features.erase(pimpl->feature_name_estimate);
      } else {
	const double score = result.prob_ + result.bound_;
	
	if (score != 0.0)
	  features[pimpl->feature_name] = score;
	else
	  features.erase(pimpl->feature_name);
      }

      if (result.oov_)
	features[pimpl->feature_name_oov] = - result.oov_;
      else
	features.erase(pimpl->feature_name_oov);      
    }

    // temporarily assigned feature function...
    
    void NGramNN::apply_predict(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				const bool final) const
    {
      // add <s>
      if (final) {
	impl_type::result_type result;
	
	pimpl->ngram_predict_score(state, result);
      }
    }
    
    void NGramNN::apply_scan(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     const int dot,
			     feature_set_type& features,
			     const bool final) const
    {
      impl_type::result_type result;
      
      pimpl->ngram_scan_score(state, edge, dot, result);
      
      if (result.prob_ != 0.0)
	features[pimpl->feature_name] = result.prob_;
      else
	features.erase(pimpl->feature_name);
      
      if (result.oov_)
	features[pimpl->feature_name_oov] = - result.oov_;
      else
	features.erase(pimpl->feature_name_oov);
    }
    
    void NGramNN::apply_complete(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 const bool final) const
    {
      // if final, add scoring for </s>
      
      if (final) {
	impl_type::result_type result;
	
	pimpl->ngram_complete_score(state, result);
	
	if (result.prob_ != 0.0)
	  features[pimpl->feature_name] = result.prob_;
	else
	  features.erase(pimpl->feature_name);
      }
    }
  };
};
