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
      
      typedef cicada::NGram      ngram_type;
      typedef cicada::NGramCache ngram_cache_type;

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
      
      typedef std::vector<double, std::allocator<double> > decay_set_type;

      typedef utils::hashmurmur<size_t> hasher_type;


      struct CacheContext
      {
	typedef utils::simple_vector<symbol_type, std::allocator<symbol_type> > phrase_type;

	phrase_type ngram;
	double logprob;
	int pos;
	
	CacheContext() : ngram(), logprob(0.0), pos(0) {}
      };
      
      struct CacheNGram
      {
	typedef utils::simple_vector<symbol_type, std::allocator<symbol_type> > phrase_type;
	
	phrase_type ngram;
	double logprob;
	
	CacheNGram() : ngram(), logprob(0.0) {}
      };
      
      typedef CacheContext cache_context_type;
      typedef CacheNGram   cache_ngram_type;
      
      typedef utils::array_power2<cache_context_type,  1024 * 128, std::allocator<cache_context_type> >  cache_context_set_type;
      typedef utils::array_power2<cache_ngram_type,    1024 * 64,  std::allocator<cache_ngram_type> >    cache_ngram_set_type;
      
    public:
      typedef boost::filesystem::path path_type;
      
      NGramImpl(const path_type& __path, const int __order)
	: ngram(__path), order(__order), cluster(0), coarse(false), no_bos_eos(false), skip_sgml_tag(false)
      {
	order = utils::bithack::min(order, ngram.index.order());
	
	initialize_decay();
	initialize_cache();
	
	id_oov = ngram.index.vocab()[vocab_type::UNK];
      }

      NGramImpl(const NGramImpl& x)
	: ngram(x.ngram),
	  order(x.order),
	  cluster(x.cluster ? &cluster_type::create(x.cluster->path()) : 0),
	  coarse(x.coarse),
	  no_bos_eos(x.no_bos_eos),
	  skip_sgml_tag(x.skip_sgml_tag),
	  feature_name(x.feature_name),
	  feature_name_oov(x.feature_name_oov),
	  id_oov(x.id_oov)
	  
      {
	initialize_decay();
	initialize_cache();
      }

      NGramImpl& operator=(const NGramImpl& x)
      {
	ngram = x.ngram;
	order = x.order;
	cluster = (x.cluster ? &cluster_type::create(x.cluster->path()) : 0);
	coarse = x.coarse;
	no_bos_eos = x.no_bos_eos;
	skip_sgml_tag = x.skip_sgml_tag;
	
	feature_name     = x.feature_name;
	feature_name_oov = x.feature_name_oov;
	id_oov           = x.id_oov;
		
	initialize_decay();
	initialize_cache();
	
	return *this;
      }

      void initialize_decay()
      {
	decays.clear();
	decays.resize(order + 1, 1.0);
	
	for (int n = 1; n < order - 1; ++ n) {
	  // 0.8 ^ (order - 1 - n)
	  for (int i = 0; i < order - 1 - n; ++ i)
	    decays[n] *= 0.8;
	}
      }

      void initialize_cache()
      {
	cache_logprob.clear();
	cache_estimate.clear();
	
	ngram_cache_logprob  = ngram_cache_type(ngram.index.order());
	ngram_cache_logbound = ngram_cache_type(ngram.index.order());
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

      template <typename Iterator1, typename Iterator2>
      inline
      bool equal_phrase(Iterator1 first1, Iterator1 last1, Iterator2 first2, Iterator2 last2) const
      {
	return std::distance(first1, last1) == std::distance(first2, last2) && std::equal(first1, last1, first2);
      }
      
      
      template <typename Iterator>
      double ngram_logprob(Iterator first, Iterator last) const
      {
	if (first == last) return 0.0;

	first = std::max(first, last - ngram.index.order());
	
	if (std::distance(first, last) <= 2)
	  return ngram.logprob(first, last);
	
	const size_type cache_pos = ngram_cache_logprob(first, last);
	if (! ngram_cache_logprob.equal_to(cache_pos, first, last)) {
	  const_cast<ngram_cache_type&>(ngram_cache_logprob).score(cache_pos) = ngram.logprob(first, last);
	  const_cast<ngram_cache_type&>(ngram_cache_logprob).assign(cache_pos, first, last);
	}
	
	return ngram_cache_logprob.score(cache_pos);
      }
      
      template <typename Iterator>
      double ngram_logbound(Iterator first, Iterator last) const
      {
	if (first == last) return 0.0;
	
	first = std::max(first, last - ngram.index.order());
	
	if (std::distance(first, last) <= 2)
	  return ngram.logbound(first, last);
	
	const size_type cache_pos = ngram_cache_logbound(first, last);
	if (! ngram_cache_logbound.equal_to(cache_pos, first, last)) {
	  const_cast<ngram_cache_type&>(ngram_cache_logbound).score(cache_pos) = ngram.logbound(first, last);
	  const_cast<ngram_cache_type&>(ngram_cache_logbound).assign(cache_pos, first, last);
	}
	
	return ngram_cache_logbound.score(cache_pos);
      }
      
      template <typename Iterator>
      double ngram_score(Iterator first, Iterator iter, Iterator last) const
      {
	if (iter == last) return 0.0;
	
	first = std::max(first, iter - order + 1);
	
	if (std::distance(iter, last) == 1) {
	  buffer_id_type& buffer_id = const_cast<buffer_id_type&>(buffer_id_impl);
	  buffer_id.clear();
	  
	  for (/**/; first != last; ++ first)
	    buffer_id.push_back(ngram.index.vocab()[*first]);
	  
	  return ngram_logprob(buffer_id.begin(), buffer_id.end());
	}
	
	const size_t cache_pos = hash_phrase(first, last, last - iter) & (cache_logprob.size() - 1);
	cache_context_type& cache = const_cast<cache_context_type&>(cache_logprob[cache_pos]);
	
	cache_context_type::phrase_type::const_iterator citer_begin = cache.ngram.begin();
	cache_context_type::phrase_type::const_iterator citer       = citer_begin + cache.pos;
	cache_context_type::phrase_type::const_iterator citer_end   = cache.ngram.end();

	if (citer_begin == citer_end || ! equal_phrase(first, iter, citer_begin, citer) || ! equal_phrase(iter, last, citer, citer_end)) {
	  cache.ngram.assign(first, last);
	  cache.pos = iter - first;
	  cache.logprob = 0.0;
	  
	  buffer_id_type& buffer_id = const_cast<buffer_id_type&>(buffer_id_impl);
	  buffer_id.clear();
	  
	  for (/**/; first != iter; ++ first)
	    buffer_id.push_back(ngram.index.vocab()[*first]);
	  
	  if (coarse) {
	    for (/**/; iter != last; ++ iter) {
	      buffer_id.push_back(ngram.index.vocab()[*iter]);
	      //cache.logprob += ngram.logbound(std::max(buffer_id.begin(), buffer_id.end() - order), buffer_id.end());
	      cache.logprob += ngram_logbound(std::max(buffer_id.begin(), buffer_id.end() - order), buffer_id.end());
	    }
	  } else {
	    for (/**/; iter != last; ++ iter) {
	      buffer_id.push_back(ngram.index.vocab()[*iter]);
	      //cache.logprob += ngram.logprob(std::max(buffer_id.begin(), buffer_id.end() - order), buffer_id.end());
	      cache.logprob += ngram_logprob(std::max(buffer_id.begin(), buffer_id.end() - order), buffer_id.end());
	    }
	  }
	}
	
	return cache.logprob;
      }
      
      template <typename Iterator>
      double ngram_estimate(Iterator first, Iterator last) const
      {
	if (first == last) return 0.0;
	if (std::distance(first, last) == 1)
	  return (vocab_type::BOS == *first ? 0.0 : ngram.logbound(first, last));

	const size_t cache_pos = hash_phrase(first, last) & (cache_estimate.size() - 1);
	cache_ngram_type& cache = const_cast<cache_ngram_type&>(cache_estimate[cache_pos]);
	
	if (! equal_phrase(first, last, cache.ngram)) {
	  cache.ngram.assign(first, last);
	  cache.logprob = 0.0;
	  
	  buffer_id_type& buffer_id = const_cast<buffer_id_type&>(buffer_id_impl);
	  buffer_id.clear();
	  
	  for (/**/; first != last; ++ first) {
	    buffer_id.push_back(ngram.index.vocab()[*first]);
	    
	    if (buffer_id.size() == 1 && vocab_type::BOS == *first) continue;
	    
	    cache.logprob += ngram_logbound(buffer_id.begin(), buffer_id.end());
	  }
	}
	  
	return cache.logprob;
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
	buffer.reserve((titer_end - titer_begin) + (order * 2) * states.size());
	
	if (states.empty()) {
	  symbol_type* context = reinterpret_cast<symbol_type*>(state);
	  std::fill(context, context + order * 2, vocab_type::EMPTY);
	  
	  // we will copy to buffer...
	  for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	    if (! skipper(*titer)) {
	      buffer.push_back(extract(*titer));
	      oov += (ngram.index.vocab()[buffer.back()] == id_oov);
	    }
	  
	  if (static_cast<int>(buffer.size()) <= context_size) {
	    std::copy(buffer.begin(), buffer.end(), context);
	    return 0.0;
	  } else {
	    buffer_type::const_iterator biter_begin = buffer.begin();
	    buffer_type::const_iterator biter_end   = buffer.end();
	    
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix = ngram.ngram_prefix(biter_begin, biter_begin + context_size);
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> suffix = ngram.ngram_suffix(biter_end - context_size, biter_end);
	    
	    std::copy(prefix.first, prefix.second, context);
	    context[prefix.second - prefix.first] = vocab_type::STAR;
	    std::copy(suffix.first, suffix.second, context + (prefix.second - prefix.first) + 1);
	    
	    return ngram_score(prefix.first, prefix.second, suffix.second);
	  }
	}

	double score = 0.0;
	
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
	    
	    const symbol_type* context = reinterpret_cast<const symbol_type*>(states[antecedent_index]);
	    const symbol_type* context_end  = std::find(context, context + order * 2, vocab_type::EMPTY);
	    const symbol_type* context_star = std::find(context, context_end, vocab_type::STAR);
	    
	    // subtract estimated score
	    //score -= ngram_estimate(context, context_star);
	    
	    buffer.insert(buffer.end(), context, context_star);
	    if (biter - biter_first >= context_size || star_first >= 0)
	      score += ngram_score(biter_first, biter, buffer.end());
	    else
	      score += ngram_score(biter_first, std::min(biter_first + context_size, buffer.end()), buffer.end());
	    biter = buffer.end();
	    
	    if (context_star != context_end) {
	      star_first = utils::bithack::branch(star_first < 0, static_cast<int>(buffer.size()) + 1, star_first);
	      star_last  = buffer.size() + 1;
	      
	      biter_first = buffer.end() + 1;
	      buffer.insert(buffer.end(), context_star, context_end);
	      biter = buffer.end();
	    }
	    
	  } else if (! skipper(*titer)) {
	    buffer.push_back(extract(*titer));
	    oov += (ngram.index.vocab()[buffer.back()] == id_oov);
	  }
	}
	
	if (biter != buffer.end()) {
	  if (biter - biter_first >= context_size || star_first >= 0)
	    score += ngram_score(biter_first, biter, buffer.end());
	  else
	    score += ngram_score(biter_first, std::min(biter_first + context_size, buffer.end()), buffer.end());
	}
	
	// construct state vector..
	symbol_type* context = reinterpret_cast<symbol_type*>(state);
	std::fill(context, context + order * 2, vocab_type::EMPTY);
	
	if (star_first >= 0) {
	  const int prefix_size = utils::bithack::min(star_first, context_size);
	  const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
	  
	  buffer_type::const_iterator biter_begin = buffer.begin();
	  buffer_type::const_iterator biter_end   = buffer.end();
	  
	  std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix = ngram.ngram_prefix(biter_begin, biter_begin + prefix_size);
	  std::pair<buffer_type::const_iterator, buffer_type::const_iterator> suffix = ngram.ngram_suffix(biter_end - suffix_size, biter_end);
	  
	  std::copy(prefix.first, prefix.second, context);
	  context[prefix.second - prefix.first] = vocab_type::STAR;
	  std::copy(suffix.first, suffix.second, context + (prefix.second - prefix.first) + 1);
	  
	  // add score from prefix.second to biter_begin + context_size
	  score += ngram_score(prefix.first, prefix.second, biter_begin + prefix_size);
	} else {
	  if (static_cast<int>(buffer.size()) <= context_size)
	    std::copy(buffer.begin(), buffer.end(), context);
	  else {
	    buffer_type::const_iterator biter_begin = buffer.begin();
	    buffer_type::const_iterator biter_end   = buffer.end();
	    
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix = ngram.ngram_prefix(biter_begin, biter_begin + context_size);
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> suffix = ngram.ngram_suffix(biter_end - context_size, biter_end);
	    
	    std::copy(prefix.first, prefix.second, context);
	    context[prefix.second - prefix.first] = vocab_type::STAR;
	    std::copy(suffix.first, suffix.second, context + (prefix.second - prefix.first) + 1);
	    
	    // add score from prefix.second to biter_begin + context_size
	    score += ngram_score(prefix.first, prefix.second, biter_begin + context_size);
	  }
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
	buffer.reserve(titer_end - titer_begin);
	
	double score = 0.0;
	for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	  if (titer->is_non_terminal()) {
	    if (! buffer.empty()) {
	      buffer_type::iterator biter_begin = buffer.begin();
	      buffer_type::iterator biter_end   = buffer.end();
	      buffer_type::iterator biter = std::min(biter_begin + context_size, biter_end);
	      
	      score += ngram_estimate(biter_begin, biter);
	      score += ngram_score(biter_begin, biter, biter_end);
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
	  buffer_type::iterator biter = std::min(biter_begin + context_size, biter_end);
	  
	  score += ngram_estimate(biter_begin, biter);
	  score += ngram_score(biter_begin, biter, biter_end);
	}

	return score;
      }

      
      double ngram_final_score(const state_ptr_type& state) const
      {
	const symbol_type* context      = reinterpret_cast<const symbol_type*>(state);
	const symbol_type* context_end  = std::find(context, context + order * 2, vocab_type::EMPTY);
	const symbol_type* context_star = std::find(context, context_end, vocab_type::STAR);
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	
	if (no_bos_eos) {
	  buffer.insert(buffer.end(), context, context_star);

	  return (! buffer.empty()
		  ? ngram_score(buffer.begin(), buffer.begin() + (buffer.front() == vocab_type::BOS), buffer.end())
		  : 0.0);
	} else {
	  buffer.push_back(vocab_type::BOS);
	  buffer.insert(buffer.end(), context, context_star);
	  
	  double score = ngram_score(buffer.begin(), buffer.begin() + 1, buffer.end());
	  
	  if (context_star != context_end) {
	    buffer.clear();
	    buffer.insert(buffer.end(), context_star + 1, context_end);
	  }
	  buffer.push_back(vocab_type::EOS);
	  
	  score += ngram_score(buffer.begin(), buffer.end() - 1, buffer.end());
	 
	  return score;
	}
      }
      
      double ngram_estimate(const state_ptr_type& state) const
      {
	const symbol_type* context = reinterpret_cast<const symbol_type*>(state);
	const symbol_type* context_end = context + order - 1;
	
	if (*context == vocab_type::EMPTY)
	  return 0.0;
	
	const symbol_type* citer_end = context;
	
	for (/**/; citer_end != context_end && *citer_end != vocab_type::STAR && *citer_end != vocab_type::EMPTY; ++ citer_end) {}
	
	return ngram_estimate(context, citer_end);
      }
      
      double ngram_predict_score(const state_ptr_type& state)
      {
	symbol_type* context = reinterpret_cast<symbol_type*>(state);
	
	std::fill(context, context + order * 2, vocab_type::EMPTY);
	
	context[0] = vocab_type::BOS;
	
	return 0.0;
      }

      double ngram_scan_score(state_ptr_type& state,
			      const edge_type& edge,
			      const int dot)
      {
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
      }
      
      double ngram_complete_score(state_ptr_type& state)
      {
	symbol_type* context = reinterpret_cast<symbol_type*>(state);
	symbol_type* context_end  = std::find(context, context + order, vocab_type::EMPTY);
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	buffer.insert(buffer.end(), context, context_end);
	buffer.push_back(vocab_type::EOS);
	
	return ngram_score(buffer.begin(), buffer.end() - 1, buffer.end());
      }

      
      decay_set_type decays;
      
      // caching...
      cache_context_set_type  cache_logprob;
      cache_ngram_set_type    cache_estimate;
      
      ngram_cache_type  ngram_cache_logprob;
      ngram_cache_type  ngram_cache_logbound;
      
      // actual buffers...
      buffer_type    buffer_impl;
      buffer_id_type buffer_id_impl;
      
      // ngrams
      ngram_type     ngram;
      int            order;

      // cluster...
      cluster_type* cluster;
      
      bool coarse;
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

      ngram_impl->no_bos_eos = no_bos_eos;
      ngram_impl->skip_sgml_tag = skip_sgml_tag;
      
      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file: " + cluster_path.string());
	
	ngram_impl->cluster = &cicada::Cluster::create(cluster_path);
      }
      
      // two contexts (order - 1) for each edge, with two separator..
      base_type::__state_size = sizeof(symbol_type) * ngram_impl->order * 2;
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
		      feature_set_type& estimates,
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
      
      if (! final) {
	const double estimate = pimpl->ngram_estimate(state);
	if (estimate != 0.0)
	  estimates[pimpl->feature_name] = estimate;
	else
	  estimates.erase(pimpl->feature_name);
      }
    }

    void NGram::apply_coarse(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features,
			     feature_set_type& estimates,
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
	
	if (! final) {
	  const double estimate = pimpl_coarse->ngram_estimate(state);
	  if (estimate != 0.0)
	    estimates[pimpl->feature_name] = estimate;
	  else
	    estimates.erase(pimpl->feature_name);
	}
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
			      feature_set_type& estimates,
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
			   feature_set_type& estimates,
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
			       feature_set_type& estimates,
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
