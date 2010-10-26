#include <stdexcept>
#include <memory>

#include "cicada/ngram.hpp"
#include "cicada/feature/ngram.hpp"
#include "cicada/parameter.hpp"
#include "cicada/symbol_vector.hpp"
#include "cicada/cluster.hpp"

#include "utils/array_power2.hpp"
#include "utils/hashmurmur.hpp"


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
      
      typedef cicada::NGram ngram_type;

      typedef cicada::Cluster cluster_type;
      
      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;
      typedef std::vector<symbol_type::id_type, std::allocator<symbol_type::id_type> > buffer_id_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;
      
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;
      
      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;

      typedef std::vector<double, std::allocator<double> > decay_set_type;

      typedef utils::hashmurmur<size_t> hasher_type;


      struct CacheContext
      {
	typedef utils::simple_vector<symbol_type, std::allocator<symbol_type> > phrase_type;

	phrase_type context;
	phrase_type ngram;
	double logprob;
	
	CacheContext() : context(), ngram(), logprob(0.0) {}
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
	: ngram(__path), order(__order), cluster(0), yield_source(false), coarse(false)
      {
	order = utils::bithack::min(order, ngram.index.order());
	
	cache_logprob.clear();
	cache_estimate.clear();
	
	initialize_decay();
      }

      NGramImpl(const NGramImpl& x)
	: ngram(x.ngram), order(x.order), cluster(x.cluster), yield_source(x.yield_source), coarse(x.coarse)
      {
	cache_logprob.clear();
	cache_estimate.clear();
	
	initialize_decay();
      }

      NGramImpl& operator=(const NGramImpl& x)
      {
	ngram = x.ngram;
	order = x.order;
	cluster = x.cluster;
	yield_source = x.yield_source;
	coarse = x.coarse;
	
	cache_logprob.clear();
	cache_estimate.clear();
	
	initialize_decay();
	
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
	  cache.logprob = 0.0;
	  
	  buffer_id_type& buffer_id = const_cast<buffer_id_type&>(buffer_id_impl);
	  buffer_id.clear();
	  buffer_id.reserve(last - first);
	  
	  for (/**/; first != iter; ++ first)
	    buffer_id.push_back(ngram.index.vocab()[*first]);
	  
	  for (/**/; iter != last; ++ iter) {
	    buffer_id.push_back(ngram.index.vocab()[*iter]);
	    
	    if (coarse)
	      cache.logprob += ngram.logbound(std::max(buffer_id.begin(), buffer_id.end() - order), buffer_id.end());
	    else
	      cache.logprob += ngram.logprob(std::max(buffer_id.begin(), buffer_id.end() - order), buffer_id.end());
	  }
	}
	
	return cache.logprob;
      }
      
      template <typename Iterator>
      double ngram_estimate(Iterator first, Iterator last) const
      {
	const size_t cache_pos = hash_phrase(first, last) & (cache_estimate.size() - 1);
	cache_ngram_type& cache = const_cast<cache_ngram_type&>(cache_estimate[cache_pos]);
	  
	if (! equal_phrase(first, last, cache.ngram)) {
	  cache.ngram.assign(first, last);
	  cache.logprob = 0.0;
	    
	  buffer_id_type& buffer_id = const_cast<buffer_id_type&>(buffer_id_impl);
	  buffer_id.clear();
	    
	  for (/**/; first != last; ++ first) {
	    buffer_id.push_back(ngram.index.vocab()[*first]);
	      
	    bool estimated = false;
	    double logbound = ngram.logbound(buffer_id.begin(), buffer_id.end(), estimated);
	      
	    if (! ngram.bound_exact && estimated && logbound < 0.0)
	      logbound *= decays[buffer_id.size()];
	      
	    cache.logprob += logbound;
	  }
	}
	  
	return cache.logprob;
      }
      
      double ngram_score(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge) const
      {
	const int context_size = order - 1;
	const rule_type& rule = *(edge.rule);
	const phrase_type& target = (yield_source ? rule.source : rule.target);
	
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
	  
	  if (cluster) {
	    for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	      if (*titer != vocab_type::EPSILON)
		buffer.push_back(cluster->operator[](*titer));
	  } else {
	    for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	      if (*titer != vocab_type::EPSILON)
		buffer.push_back(*titer);
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
	
	phrase_span_set_type& phrase_spans = const_cast<phrase_span_set_type&>(phrase_spans_impl);
	phrase_spans.clear();
	target.terminals(std::back_inserter(phrase_spans));
	
	double score = 0.0;
	
	int star_first = -1;
	int star_last  = -1;
	
	if (phrase_spans.front().second - phrase_spans.front().first >= order)
	  score += ngram_score(phrase_spans.front().first, phrase_spans.front().first + context_size, phrase_spans.front().second);
	
	buffer.insert(buffer.end(), phrase_spans.front().first, phrase_spans.front().second);
	
	buffer_type::const_iterator biter_first = buffer.begin();
	
	phrase_span_set_type::const_iterator siter_begin = phrase_spans.begin();
	phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	for (phrase_span_set_type::const_iterator siter = siter_begin + 1; siter != siter_end; ++ siter) {
	  const phrase_span_type& span = *siter;
	  
	  int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	  if (antecedent_index < 0)
	    antecedent_index = siter - (siter_begin + 1);
	  
	  const symbol_type* context = reinterpret_cast<const symbol_type*>(states[antecedent_index]);
	  const symbol_type* context_end  = std::find(context, context + order * 2, vocab_type::EMPTY);
	  const symbol_type* context_star = std::find(context, context_end, vocab_type::STAR);
	  
	  buffer_type::const_iterator biter = buffer.end();
	  
	  buffer.insert(buffer.end(), context, context_star);
	  
	  buffer_type::const_iterator biter_end = buffer.end();

	  if (biter - biter_first >= context_size || star_first >= 0)
	    score += ngram_score(biter_first, biter, biter_end);
	  else
	    score += ngram_score(biter_first, std::min(biter_first + context_size, biter_end), biter_end);
	  
	  // insert star!
	  if (context_star != context_end) {
	    biter_first = buffer.end() + 1;

	    if (star_first < 0)
	      star_first = buffer.size() + 1;
	    star_last = buffer.size() + 1;
	    
	    buffer.insert(buffer.end(), context_star, context_end);
	  }
	  
	  // use of this edge's terminals
	  {
	    buffer_type::const_iterator biter = buffer.end();
	    
	    if (cluster) {
	      for (phrase_type::const_iterator titer = span.first; titer != span.second; ++ titer)
		buffer.push_back(cluster->operator[](*titer));
	    } else
	      buffer.insert(buffer.end(), span.first, span.second);
	    
	    buffer_type::const_iterator biter_end = buffer.end();
	    
	    if (biter - biter_first >= context_size || star_first >= 0)
	      score += ngram_score(biter_first, biter, biter_end);
	    else
	      score += ngram_score(biter_first, std::min(biter_first + context_size, biter_end), biter_end);
	  }
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

      double ngram_estimate(const edge_type& edge) const
      {
	const int context_size = order - 1;
	const rule_type& rule = *(edge.rule);
	const phrase_type& target = (yield_source ? rule.source : rule.target);
	
	phrase_type::const_iterator titer_begin = target.begin();
	phrase_type::const_iterator titer_end   = target.end();
	
	// we will reserve enough space so that buffer's memory will not be re-allocated.
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	buffer.reserve(titer_end - titer_begin);
	
	double score = 0.0;
	if (cluster) {
	  for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) 
	    if (titer->is_non_terminal()) {
	      if (! buffer.empty()) {
		buffer_type::iterator biter_begin = buffer.begin();
		buffer_type::iterator biter_end   = buffer.end();
		buffer_type::iterator biter = std::min(biter_begin + context_size, biter_end);
		
		score += ngram_estimate(biter_begin, biter);
		score += ngram_score(biter_begin, biter, biter_end);
	      }
	      
	      buffer.clear();
	    } else if (*titer != vocab_type::EPSILON)
	      buffer.push_back(cluster->operator[](*titer));
	} else {
	  for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) 
	    if (titer->is_non_terminal()) {
	      if (! buffer.empty()) {
		buffer_type::iterator biter_begin = buffer.begin();
		buffer_type::iterator biter_end   = buffer.end();
		buffer_type::iterator biter = std::min(biter_begin + context_size, biter_end);
		
		score += ngram_estimate(biter_begin, biter);
		score += ngram_score(biter_begin, biter, biter_end);
	      }
	      
	      buffer.clear();
	    } else if (*titer != vocab_type::EPSILON)
	      buffer.push_back(*titer);
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
      
      double ngram_estimate(const state_ptr_type& state) const
      {
	const symbol_type* context = reinterpret_cast<const symbol_type*>(state);
	const symbol_type* context_end = context + order - 1;
	
	if (*context == vocab_type::EMPTY)
	  return 0.0;
	
	const symbol_type* citer_end = context;
	
	for (/**/; citer_end != context_end && *citer_end != vocab_type::STAR && *citer_end != vocab_type::EMPTY; ++ citer_end);
	
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
	const phrase_type& phrase = (yield_source ? rule.source : rule.target);
	
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
      
      phrase_span_set_type phrase_spans_impl;
      
      // actual buffers...
      buffer_type    buffer_impl;
      buffer_id_type buffer_id_impl;
      
      // ngrams
      ngram_type     ngram;
      int            order;

      // cluster...
      cluster_type* cluster;
      
      // yield/coarse
      bool yield_source;
      bool coarse;
    };
    
    
    NGram::NGram(const std::string& parameter)
      : pimpl(0), pimpl_coarse(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "ngram")
	throw std::runtime_error("is this really ngram feature function? " + parameter);

      path_type   path;
      int         order = 3;
      path_type   cluster_path;
      
      path_type   coarse_path;
      int         coarse_order = 0;
      path_type   coarse_cluster_path;
      
      bool        yield_source = false;
      bool        yield_target = false;
      std::string name;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "file") == 0)
	  path = piter->second;
	else if (strcasecmp(piter->first.c_str(), "cluster") == 0)
	  cluster_path = piter->second;
	else if (strcasecmp(piter->first.c_str(), "order") == 0)
	  order = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "coarse-file") == 0)
	  coarse_path = piter->second;
	else if (strcasecmp(piter->first.c_str(), "coarse-order") == 0)
	  coarse_order = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "coarse-cluster") == 0)
	  coarse_cluster_path = piter->second;
	else if (strcasecmp(piter->first.c_str(), "name") == 0)
	  name = piter->second;
	else if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	  const std::string& yield = piter->second;
	  
	  if (strcasecmp(yield.c_str(), "source") == 0)
	    yield_source = true;
	  else if (strcasecmp(yield.c_str(), "target") == 0)
	    yield_target = true;
	  else
	    throw std::runtime_error("unknown parameter: " + parameter);
	} else
	  std::cerr << "WARNING: unsupported parameter for ngram: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (path.empty())
	throw std::runtime_error("no ngram file? " + path.file_string());
      
      if (order <= 0)
	throw std::runtime_error("invalid ngram order: " + boost::lexical_cast<std::string>(order));

      if (yield_source && yield_target)
	throw std::runtime_error("you cannot specify both source/target yield");
      
      if (coarse_order > order)
	throw std::runtime_error("invalid coarse order: coarse-order <= order");
      if (! coarse_path.empty() && ! boost::filesystem::exists(coarse_path))
	throw std::runtime_error("no coarse ngram language model? " + coarse_path.file_string());
      
      
      std::auto_ptr<impl_type> ngram_impl(new impl_type(path, order));
      
      // set up yield..
      ngram_impl->yield_source = yield_source;
      
      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file: " + cluster_path.file_string());
	
	ngram_impl->cluster = &cicada::Cluster::create(cluster_path);
      }

      
      // two contexts (order - 1) for each edge, with two separator..
      base_type::__state_size = sizeof(symbol_type) * ngram_impl->order * 2;
      base_type::__feature_name = (name.empty() ? std::string("ngram") : name);
      
      pimpl = ngram_impl.release();

      // ...
      if (coarse_order > 0 || ! coarse_path.empty()) {
	
	if (coarse_order <= 0)
	  throw std::runtime_error("coarse order must be non-zero!");
	
	if (! coarse_path.empty()) {
	  std::auto_ptr<impl_type> ngram_impl(new impl_type(coarse_path, coarse_order));
	  
	  // set up yield..
	  ngram_impl->yield_source = yield_source;
	  
	  if (! coarse_cluster_path.empty()) {
	    if (! boost::filesystem::exists(coarse_cluster_path))
	      throw std::runtime_error("no cluster file: " + coarse_cluster_path.file_string());
	    
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
      double score = pimpl->ngram_score(state, states, edge);
      if (final)
	score += pimpl->ngram_final_score(state);
      
      if (score != 0.0)
	features[base_type::feature_name()] = score;
      else
	features.erase(base_type::feature_name());
      
      if (! final) {
	const double estimate = pimpl->ngram_estimate(state);
	if (estimate != 0.0)
	  estimates[base_type::feature_name()] = estimate;
	else
	  estimates.erase(base_type::feature_name());
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
	double score = pimpl_coarse->ngram_score(state, states, edge);
	if (final)
	  score += pimpl_coarse->ngram_final_score(state);
      
	if (score != 0.0)
	  features[base_type::feature_name()] = score;
	else
	  features.erase(base_type::feature_name());
      
	if (! final) {
	  const double estimate = pimpl_coarse->ngram_estimate(state);
	  if (estimate != 0.0)
	    estimates[base_type::feature_name()] = estimate;
	  else
	    estimates.erase(base_type::feature_name());
	}
      } else {
	// state-less...
	const double score = pimpl->ngram_estimate(edge);
	
	if (score != 0.0)
	  features[base_type::feature_name()] = score;
	else
	  features.erase(base_type::feature_name());
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
	features[base_type::feature_name()] = score;
      else
	features.erase(base_type::feature_name());
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
	  features[base_type::feature_name()] = score;
	else
	  features.erase(base_type::feature_name());
      }
    }


  };
};
