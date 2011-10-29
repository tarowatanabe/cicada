//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <memory>

#include "cicada/ngram.hpp"
#include "cicada/ngram_cache.hpp"
#include "cicada/feature/ngram.hpp"
#include "cicada/parameter.hpp"
#include "cicada/symbol_vector.hpp"
#include "cicada/cluster.hpp"

#include "utils/vector2.hpp"
#include "utils/array_power2.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"

namespace cicada
{
  namespace feature
  {
    class NGramImpl : public utils::hashmurmur<size_t>
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Vocab  vocab_type;
      
      typedef cicada::NGram                           ngram_type;

      typedef ngram_type::state_type ngram_state_type;
      typedef std::pair<ngram_state_type, double> state_score_type;
      
      typedef cicada::NGramCache<symbol_type, state_score_type> ngram_cache_type;
      
      typedef cicada::Cluster cluster_type;
      
      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;
      typedef std::vector<symbol_type::id_type, std::allocator<symbol_type::id_type> > buffer_id_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_set_type::feature_type feature_type;
      
      typedef rule_type::symbol_set_type phrase_type;
      
      typedef utils::hashmurmur<size_t> hasher_type;


      struct CacheContext
      {
	typedef utils::simple_vector<symbol_type, std::allocator<symbol_type> > phrase_type;

	ngram_state_type state;
	phrase_type      ngram;
	state_score_type score;
	
	CacheContext() : state(), ngram(), score() {}
      };
      
      typedef CacheContext cache_context_type;
      typedef utils::array_power2<cache_context_type,  1024 * 128, std::allocator<cache_context_type> >  cache_context_set_type;
      
    public:
      typedef boost::filesystem::path path_type;
      
      NGramImpl(const path_type& __path, const int __order)
	: ngram(__path), order(__order), cluster(0), coarse(false), approximate(false), no_bos_eos(false), skip_sgml_tag(false)
      {
	order = utils::bithack::min(order, ngram.index.order());
	
	initialize_cache();
	
	id_oov = ngram.index.vocab()[vocab_type::UNK];
      }

      NGramImpl(const NGramImpl& x)
	: ngram(x.ngram),
	  order(x.order),
	  cluster(x.cluster ? &cluster_type::create(x.cluster->path()) : 0),
	  coarse(x.coarse),
	  approximate(x.approximate),
	  no_bos_eos(x.no_bos_eos),
	  skip_sgml_tag(x.skip_sgml_tag),
	  feature_name(x.feature_name),
	  feature_name_oov(x.feature_name_oov),
	  id_oov(x.id_oov)
	  
      {
	initialize_cache();
      }

      NGramImpl& operator=(const NGramImpl& x)
      {
	ngram = x.ngram;
	order = x.order;
	cluster = (x.cluster ? &cluster_type::create(x.cluster->path()) : 0);
	coarse = x.coarse;
	approximate = x.approximate;
	no_bos_eos = x.no_bos_eos;
	skip_sgml_tag = x.skip_sgml_tag;
	
	feature_name     = x.feature_name;
	feature_name_oov = x.feature_name_oov;
	id_oov           = x.id_oov;
		
	initialize_cache();
	
	return *this;
      }
      
      void initialize_cache()
      {
	cache_logprob.clear();
	cache_estimate = ngram_cache_type(ngram.index.order());
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
      
      struct __ngram_score_logprob
      {
	const ngram_type& ngram;

	__ngram_score_logprob(const ngram_type& __ngram) : ngram(__ngram) {}
	
	template <typename Word>
	state_score_type operator()(const ngram_state_type& state, const Word& word, const bool backoffed, const int max_order) const
	{
	  return ngram.logprob(state, word, backoffed, max_order);
	}
      };

      struct __ngram_score_logbound
      {
	const ngram_type& ngram;

	__ngram_score_logbound(const ngram_type& __ngram) : ngram(__ngram) {}
	
	template <typename Word>
	state_score_type operator()(const ngram_state_type& state, const Word& word, const bool backoffed, const int max_order) const
	{
	  return ngram.logbound(state, word, backoffed, max_order);
	}
      };

      template <typename Iterator>
      state_score_type ngram_score(const ngram_state_type& state, Iterator first, Iterator last) const
      {
	if (coarse)
	  return ngram_score(state, first, last, __ngram_score_logbound(ngram));
	else
	  return ngram_score(state, first, last, __ngram_score_logprob(ngram));
      }
      
      template <typename Iterator, typename Scorer>
      state_score_type ngram_score(ngram_state_type state, Iterator first, Iterator last, Scorer scorer) const
      {
	const size_type length = std::distance(first, last);
	
	if (length == 0)
	  return state_score_type(state, 0.0);
	else if (length == 1)
	  return scorer(state, *first, ngram.index.order(state) + 1 < order, order);
	
	const size_t cache_pos = hash_phrase(first, last, state.value()) & (cache_logprob.size() - 1);
	cache_context_type& cache = const_cast<cache_context_type&>(cache_logprob[cache_pos]);
	
	if (cache.state != state || ! equal_phrase(first, last, cache.ngram)) {
	  cache.state = state;
	  cache.ngram.assign(first, last);
	  
	  ngram_type::logprob_type score = 0.0;
	  for (/**/; first != last; ++ first) {
	    const state_score_type result = scorer(state, *first, ngram.index.order(state) + 1 < order, order);
	    
	    state = result.first;
	    score += result.second;
	  }
	  
	  cache.score = std::make_pair(state, score);
	}
	
	return cache.score;
      }
      
      template <typename Iterator>
      state_score_type ngram_estimate(Iterator first, Iterator last) const
      {
	if (approximate)
	  return ngram_estimate(first, last, __ngram_score_logprob(ngram));
	else
	  return ngram_estimate(first, last, __ngram_score_logbound(ngram));
      }
      
      template <typename Iterator, typename Scorer>
      state_score_type ngram_estimate(Iterator first, Iterator last, Scorer scorer) const
      {
	const size_type length = std::distance(first, last);
	
	if (length == 0)
	  return state_score_type(ngram_state_type(), 0.0);
	else if (length <= 2) {
	  int bound_order = 0;
	  state_score_type state_score(ngram_state_type(), 0.0);
	  
	  if (*first == vocab_type::BOS) {
	    state_score.first = ngram.index.next(state_score.first, vocab_type::BOS);
	    ++ first;
	    ++ bound_order;
	  }
	  
	  for (/**/; first != last; ++ first, ++ bound_order) {
	    const state_score_type result = scorer(state_score.first, *first, ngram.index.order(state_score.first) < bound_order, order);
	    
	    state_score.first   = result.first;
	    state_score.second += result.second;
	  }
	  
	  return state_score;
	}
	
	
	const size_type cache_pos = cache_estimate(first, last);
	
	if (! cache_estimate.equal_to(cache_pos, first, last)) {
	  ngram_cache_type& cache = const_cast<ngram_cache_type&>(cache_estimate);
	  cache.assign(cache_pos, first, last);
	  
	  int bound_order = 0;
	  state_score_type state_score(ngram_state_type(), 0.0);
	  
	  if (*first == vocab_type::BOS) {
	    state_score.first = ngram.index.next(state_score.first, vocab_type::BOS);
	    ++ first;
	    ++ bound_order;
	  }
	  
	  for (/**/; first != last; ++ first, ++ bound_order) {
	    const state_score_type result = scorer(state_score.first, *first, ngram.index.order(state_score.first) < bound_order, order);
	    
	    state_score.first   = result.first;
	    state_score.second += result.second;
	  }
	  
	  cache.score(cache_pos) = state_score;
	}
	
	return cache_estimate.score(cache_pos);
      }

      struct extract_cluster
      {
	extract_cluster(const cluster_type* __cluster): cluster(__cluster) {}
	
	const cluster_type* cluster;

	symbol_type operator()(const symbol_type& word) const
	{
	  return cluster->operator[](word);
	}
      };

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
      
      double ngram_score(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 int& oov) const
      {
	if (cluster) {
	  if (skip_sgml_tag)
	    return ngram_score(state, states, edge, oov, extract_cluster(cluster), skipper_sgml());
	  else
	    return ngram_score(state, states, edge, oov, extract_cluster(cluster), skipper_epsilon());
	} else {
	  if (skip_sgml_tag)
	    return ngram_score(state, states, edge, oov, extract_word(), skipper_sgml());
	  else
	    return ngram_score(state, states, edge, oov, extract_word(), skipper_epsilon());
	}
      }

      double ngram_estimate(const edge_type& edge, int& oov) const
      {
	if (cluster) {
	  if (skip_sgml_tag)
	    return ngram_estimate(edge, oov, extract_cluster(cluster), skipper_sgml());
	  else
	    return ngram_estimate(edge, oov, extract_cluster(cluster), skipper_epsilon());
	} else {
	  if (skip_sgml_tag)
	    return ngram_estimate(edge, oov, extract_word(), skipper_sgml());
	  else
	    return ngram_estimate(edge, oov, extract_word(), skipper_epsilon());
	}
      }

      
      template <typename Extract, typename Skipper>
      double ngram_score(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 int& oov,
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
	buffer.clear();
	
	static const ngram_state_type state_root    = ngram_state_type();
	static const ngram_state_type state_invalid = ngram_state_type(size_type(0), size_type(-1));

	ngram_state_type* ngram_state = reinterpret_cast<ngram_state_type*>(state);
	symbol_type*      context     = reinterpret_cast<symbol_type*>(ngram_state + 1);
	symbol_type*      context_end = context + order;
	
	if (states.empty()) {
	  // we will copy to buffer...
	  for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	    if (! skipper(*titer)) {
	      buffer.push_back(extract(*titer));
	      oov += (ngram.index.vocab()[buffer.back()] == id_oov);
	    }
	  
	  if (static_cast<int>(buffer.size()) <= context_size) {
	    const state_score_type state_bound = ngram_estimate(buffer.begin(), buffer.end());
	    
	    *ngram_state = state_invalid;
	    std::copy(buffer.begin(), buffer.end(), context);
	    std::fill(context + buffer.size(), context_end, vocab_type::EMPTY);
	    
	    return state_bound.second;
	  } else {
	    buffer_type::const_iterator biter_begin = buffer.begin();
	    buffer_type::const_iterator biter_end   = buffer.end();
	    
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix = ngram.prefix(biter_begin, biter_begin + context_size);
	    
	    const state_score_type state_bound = ngram_estimate(prefix.first, prefix.second);
	    const state_score_type state_score = ngram_score(state_bound.first, prefix.second, biter_end);
	    
	    *ngram_state = state_score.first;
	    std::copy(prefix.first, prefix.second, context);
	    std::fill(context + (prefix.second - prefix.first), context_end, vocab_type::EMPTY);
	    
	    return state_bound.second + state_score.second;
	  }
	}
	
	ngram_state_type state_rule = state_invalid;
	double score = 0.0;
	
	int non_terminal_pos = 0;
	for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	  if (titer->is_non_terminal()) {
	    const int __non_terminal_index = titer->non_terminal_index();
	    const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	    ++ non_terminal_pos;
	    
	    const ngram_state_type* ngram_state_antecedent = reinterpret_cast<const ngram_state_type*>(states[antecedent_index]);
	    const symbol_type*      context_antecedent     = reinterpret_cast<const symbol_type*>(ngram_state_antecedent + 1);
	    const symbol_type*      context_antecedent_end = std::find(context_antecedent, context_antecedent + order, vocab_type::EMPTY);
	    
	    // subtract estimated prefix score
	    score -= ngram_estimate(context_antecedent, context_antecedent_end).second;
	    
	    buffer.insert(buffer.end(), context_antecedent, context_antecedent_end);

	    if (! buffer.empty()) {
	      if (state_rule != state_invalid) {
		// state_rule is "valid". we will compute anyway
		const state_score_type state_score = ngram_score(state_rule, buffer.begin(), buffer.end());
		
		state_rule = state_score.first;
		score     += state_score.second;
		
		buffer.clear();
	      } else if (static_cast<int>(buffer.size()) > context_size || *ngram_state_antecedent != state_invalid) {
		// state_rule is invalid, but we have enough context or we have 'star' at antecedent
		
		if (static_cast<int>(buffer.size()) <= context_size) {
		  const state_score_type state_bound = ngram_estimate(buffer.begin(), buffer.end());
		  
		  std::copy(buffer.begin(), buffer.end(), context);
		  std::fill(context + buffer.size(), context_end, vocab_type::EMPTY);
		  
		  score += state_bound.second;
		  
		} else {
		  buffer_type::const_iterator biter_begin = buffer.begin();
		  buffer_type::const_iterator biter_end   = buffer.end();
		  
		  std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix = ngram.prefix(biter_begin, biter_begin + context_size);
		  
		  const state_score_type state_bound = ngram_estimate(prefix.first, prefix.second);
		  const state_score_type state_score = ngram_score(state_bound.first, prefix.second, biter_end);
		  
		  state_rule = state_score.first;
		  std::copy(prefix.first, prefix.second, context);
		  std::fill(context + (prefix.second - prefix.first), context_end, vocab_type::EMPTY);
		  
		  score += state_bound.second + state_score.second;
		}
		
		buffer.clear();
	      }
	    }
	    
	    if (*ngram_state_antecedent != state_invalid)
	      state_rule = *ngram_state_antecedent;
	    
	  } else if (! skipper(*titer)) {
	    buffer.push_back(extract(*titer));
	    oov += (ngram.index.vocab()[buffer.back()] == id_oov);
	  }
	}
	
	if (buffer.empty()) {
	  *ngram_state = state_rule;
	  
	  if (state_rule == state_invalid)
	    std::fill(context, context_end, vocab_type::EMPTY);
	} else if (state_rule != state_invalid) {
	  const state_score_type state_score = ngram_score(state_rule, buffer.begin(), buffer.end());
	  
	  *ngram_state = state_score.first;
	  
	  score += state_score.second;
	} else if (static_cast<int>(buffer.size()) <= context_size) {
	  const state_score_type state_bound = ngram_estimate(buffer.begin(), buffer.end());
	  
	  *ngram_state = state_invalid;
	  std::copy(buffer.begin(), buffer.end(), context);
	  std::fill(context + buffer.size(), context_end, vocab_type::EMPTY);

	  score += state_bound.second;
	} else {
	  buffer_type::const_iterator biter_begin = buffer.begin();
	  buffer_type::const_iterator biter_end   = buffer.end();
	  
	  std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix = ngram.prefix(biter_begin, biter_begin + context_size);
	  
	  const state_score_type state_bound = ngram_estimate(prefix.first, prefix.second);
	  const state_score_type state_score = ngram_score(state_bound.first, prefix.second, biter_end);
	  
	  *ngram_state = state_score.first;
	  std::copy(prefix.first, prefix.second, context);
	  std::fill(context + (prefix.second - prefix.first), context_end, vocab_type::EMPTY);

	  score += state_bound.second + state_score.second;
	}

	return score;
      }
      
      template <typename Extract, typename Skipper>
      double ngram_estimate(const edge_type& edge, int& oov, Extract extract, Skipper skipper) const
      {
	const int context_size = order - 1;
	const rule_type& rule = *(edge.rule);
	const phrase_type& target = rule.rhs;
	
	phrase_type::const_iterator titer_begin = target.begin();
	phrase_type::const_iterator titer_end   = target.end();
	
	// we will reserve enough space so that buffer's memory will not be re-allocated.
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	
	double score = 0.0;
	for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	  if (titer->is_non_terminal()) {
	    if (! buffer.empty()) {
	      buffer_type::iterator biter_begin = buffer.begin();
	      buffer_type::iterator biter_end   = buffer.end();
	      buffer_type::iterator biter       = std::min(biter_begin + context_size, biter_end);

	      const state_score_type state_bound = ngram_estimate(biter_begin, biter);
	      const state_score_type state_score = ngram_score(state_bound.first, biter, biter_end);
	      
	      score += state_bound.second + state_score.second;
	    }
	    
	    buffer.clear();
	  } else if (! skipper(*titer)) {
	    buffer.push_back(extract(*titer));
	    oov += (ngram.index.vocab()[buffer.back()] == id_oov);
	  }
	}
	
	if (! buffer.empty()) {
	  buffer_type::iterator biter_begin = buffer.begin();
	  buffer_type::iterator biter_end   = buffer.end();
	  buffer_type::iterator biter       = std::min(biter_begin + context_size, biter_end);
	  
	  const state_score_type state_bound = ngram_estimate(biter_begin, biter);
	  const state_score_type state_score = ngram_score(state_bound.first, biter, biter_end);
	  
	  score += state_bound.second + state_score.second;
	}

	return score;
      }

      
      double ngram_final_score(const state_ptr_type& state) const
      {
	static const ngram_state_type state_root    = ngram_state_type();
	static const ngram_state_type state_invalid = ngram_state_type(size_type(0), size_type(-1));

	const ngram_state_type* ngram_state = reinterpret_cast<const ngram_state_type*>(state);
	const symbol_type*      context     = reinterpret_cast<const symbol_type*>(ngram_state + 1);
	const symbol_type*      context_end = std::find(context, context + order, vocab_type::EMPTY);
	
	if (no_bos_eos) {
	  if (context == context_end)
	    return 0.0;
	  else {
	    const state_score_type state_bound = ngram_estimate(context, context_end);
	    const state_score_type state_score = (*context == vocab_type::BOS
						  ? ngram_score(ngram.index.next(ngram_state_type(), *context), context + 1, context_end)
						  : ngram_score(ngram_state_type(), context, context_end));
	    
	    return state_score.second - state_bound.second;
	  }
	} else {
	  const state_score_type state_bound = ngram_estimate(context, context_end);
	  const state_score_type state_score = ngram_score(ngram.index.next(ngram_state_type(), vocab_type::BOS), context, context_end);
	  
	  return (state_score.second - state_bound.second
		  + ngram_score(*ngram_state != state_invalid ? *ngram_state : state_score.first, &vocab_type::EOS, (&vocab_type::EOS) + 1).second);
	}
      }
            
      double ngram_predict_score(const state_ptr_type& state)
      {
#if 0
	symbol_type* context = reinterpret_cast<symbol_type*>(state);
	
	std::fill(context, context + order * 2, vocab_type::EMPTY);
	
	context[0] = vocab_type::BOS;
	
	return 0.0;
#endif
	return 0.0;
      }

      double ngram_scan_score(state_ptr_type& state,
			      const edge_type& edge,
			      const int dot)
      {
#if 0
	const rule_type& rule = *(edge.rule);
	const phrase_type& phrase = rule.rhs;
	
	symbol_type* context = reinterpret_cast<symbol_type*>(state);
	symbol_type* context_end  = std::find(context, context + order, vocab_type::EMPTY);
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	buffer.reserve(order + phrase.size());
	buffer.insert(buffer.end(), context, context_end);
	
	buffer_type::iterator biter = buffer.end();
	
	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin() + dot; piter != piter_end && ! piter->is_non_terminal(); ++ piter)
	  if (*piter != vocab_type::EPSILON)
	    buffer.push_back(*piter);
	
	const double score = ngram_score(buffer.begin(), biter, buffer.end());
	
	std::pair<buffer_type::iterator, buffer_type::iterator> suffix = ngram.ngram_suffix(buffer.begin(), buffer.end());
	
	std::copy(suffix.first, suffix.second, context);
	std::fill(context + (suffix.second - suffix.first), context + order * 2, vocab_type::EMPTY);
	
	return score;
#endif
	return 0.0;
      }
      
      double ngram_complete_score(state_ptr_type& state)
      {
#if 0
	symbol_type* context = reinterpret_cast<symbol_type*>(state);
	symbol_type* context_end  = std::find(context, context + order, vocab_type::EMPTY);
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	buffer.insert(buffer.end(), context, context_end);
	buffer.push_back(vocab_type::EOS);
	
	return ngram_score(buffer.begin(), buffer.end() - 1, buffer.end());
#endif
	return 0.0;
      }

      
      // caching...
      cache_context_set_type cache_logprob;
      ngram_cache_type       cache_estimate;
      
      // actual buffers...
      buffer_type    buffer_impl;
      buffer_id_type buffer_id_impl;
      
      // ngrams
      ngram_type     ngram;
      int            order;

      // cluster...
      cluster_type* cluster;
      
      bool coarse;
      bool approximate;
      bool no_bos_eos;
      bool skip_sgml_tag;
      
      // names...
      feature_type feature_name;
      feature_type feature_name_oov;
      symbol_type::id_type id_oov;
    };
    
    
    NGram::NGram(const std::string& parameter)
      : pimpl(0), pimpl_coarse(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (utils::ipiece(param.name()) != "ngram")
	throw std::runtime_error("is this really ngram feature function? " + parameter);

      path_type   path;
      int         order = 3;
      path_type   cluster_path;
      bool        approximate = false;
      bool        skip_sgml_tag = false;
      bool        no_bos_eos = false;
      
      path_type   coarse_path;
      int         coarse_order = 0;
      path_type   coarse_cluster_path;
      
      std::string name;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else if (utils::ipiece(piter->first) == "cluster")
	  cluster_path = piter->second;
	else if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "approximate")
	  approximate = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "no-bos-eos")
	  no_bos_eos = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "coarse-file")
	  coarse_path = piter->second;
	else if (utils::ipiece(piter->first) == "coarse-order")
	  coarse_order = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "coarse-cluster")
	  coarse_cluster_path = piter->second;
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
      
      std::auto_ptr<impl_type> ngram_impl(new impl_type(path, order));

      ngram_impl->approximate = approximate;
      ngram_impl->no_bos_eos = no_bos_eos;
      ngram_impl->skip_sgml_tag = skip_sgml_tag;
      
      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file: " + cluster_path.string());
	
	ngram_impl->cluster = &cicada::Cluster::create(cluster_path);
      }
      
      // two contexts (order - 1) for each edge, with two separator..
      base_type::__state_size = sizeof(impl_type::ngram_state_type) + sizeof(symbol_type) * ngram_impl->order;
      base_type::__feature_name = (name.empty() ? std::string("ngram") : name);
      
      ngram_impl->feature_name     = base_type::__feature_name;
      ngram_impl->feature_name_oov = static_cast<const std::string&>(base_type::__feature_name) + ":oov-penalty";
      
      pimpl = ngram_impl.release();

      // ...
      if (coarse_order > 0 || ! coarse_path.empty()) {
	
	if (coarse_order <= 0)
	  throw std::runtime_error("coarse order must be non-zero!");
	
	if (! coarse_path.empty()) {
	  std::auto_ptr<impl_type> ngram_impl(new impl_type(coarse_path, coarse_order));

	  ngram_impl->approximate = approximate;
	  ngram_impl->no_bos_eos = no_bos_eos;
	  ngram_impl->skip_sgml_tag = skip_sgml_tag;	  
	  
	  if (! coarse_cluster_path.empty()) {
	    if (! boost::filesystem::exists(coarse_cluster_path))
	      throw std::runtime_error("no cluster file: " + coarse_cluster_path.string());
	    
	    ngram_impl->cluster = &cicada::Cluster::create(coarse_cluster_path);
	  }
	  
	  pimpl_coarse = ngram_impl.release();
	} else {
	  std::auto_ptr<impl_type> ngram_impl(new impl_type(*pimpl));
	  ngram_impl->order = coarse_order;
	  ngram_impl->coarse = true;	  
	  
	  pimpl_coarse = ngram_impl.release();
	}
      }
    }
    
    NGram::~NGram()
    {
      std::auto_ptr<impl_type> tmp(pimpl);
      if (pimpl_coarse)
	std::auto_ptr<impl_type> tmp_coarse(pimpl_coarse);
    }
    
    NGram::NGram(const NGram& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl)),
	pimpl_coarse(x.pimpl_coarse ? new impl_type(*x.pimpl_coarse) : 0)
    {}
    
    NGram& NGram::operator=(const NGram& x)
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
    
    
    void NGram::apply(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      feature_set_type& features,
		      const bool final) const
    {
      int oov = 0;
      double score = pimpl->ngram_score(state, states, edge, oov);
      if (final)
	score += pimpl->ngram_final_score(state);
      
      if (score != 0.0)
	features[pimpl->feature_name] = score;
      else
	features.erase(pimpl->feature_name);

      if (oov)
	features[pimpl->feature_name_oov] = - oov;
      else
	features.erase(pimpl->feature_name_oov);
    }

    void NGram::apply_coarse(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features,
			     const bool final) const
    {
      if (pimpl_coarse) {
	int oov = 0;
	double score = pimpl_coarse->ngram_score(state, states, edge, oov);
	if (final)
	  score += pimpl_coarse->ngram_final_score(state);
      
	if (score != 0.0)
	  features[pimpl->feature_name] = score;
	else
	  features.erase(pimpl->feature_name);
	
	if (oov)
	  features[pimpl->feature_name_oov] = - oov;
	else
	  features.erase(pimpl->feature_name_oov);
      } else {
	// state-less.
	int oov = 0;
	const double score = pimpl->ngram_estimate(edge, oov);
	
	if (score != 0.0)
	  features[pimpl->feature_name] = score;
	else
	  features.erase(pimpl->feature_name);

	if (oov)
	  features[pimpl->feature_name_oov] = - oov;
	else
	  features.erase(pimpl->feature_name_oov);
      }
    }

    // temporarily assigned feature function...
    
    void NGram::apply_predict(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      const bool final) const
    {
      // add <s>
      if (final)
	pimpl->ngram_predict_score(state);
    }
    
    void NGram::apply_scan(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge,
			   const int dot,
			   feature_set_type& features,
			   const bool final) const
    {
      const double score = pimpl->ngram_scan_score(state, edge, dot);
      
      if (score != 0.0)
	features[pimpl->feature_name] = score;
      else
	features.erase(pimpl->feature_name);
    }
    
    void NGram::apply_complete(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       const bool final) const
    {
      // if final, add scoring for </s>
      if (final) {
	const double score = pimpl->ngram_complete_score(state);
	
	if (score != 0.0)
	  features[pimpl->feature_name] = score;
	else
	  features.erase(pimpl->feature_name);
      }
    }


  };
};
