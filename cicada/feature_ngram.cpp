#include <stdexcept>
#include <memory>

#include "cicada/ngram.hpp"
#include "cicada/feature_ngram.hpp"
#include "cicada/parameter.hpp"
#include "cicada/symbol_vector.hpp"

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
      
      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;
      typedef std::vector<symbol_type::id_type, std::allocator<symbol_type::id_type> > buffer_id_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef std::vector<double, std::allocator<double> > decay_set_type;

      typedef SymbolVector phrase_type;

      typedef utils::hashmurmur<size_t> hasher_type;
            
    public:
      typedef boost::filesystem::path path_type;
      
      NGramImpl(const path_type& __path, const int __order)
	: ngram(__path), order(__order)
      {
	order = utils::bithack::min(order, ngram.index.order());

	initialize_decay();
      }

      void initialize_decay()
      {
	decays.clear();
	decays.resize(order + 1, 1.0);
	
	for (int n = 1; n < order - 1; ++ n) {
	  // 0.7 ^ (order - 1 - n)
	  for (int i = 0; i < order - 1 - n; ++ i)
	    decays[n] *= 0.7;
	}
      }

      template <typename Iterator>
      inline
      size_t hash_phrase(Iterator first, Iterator last, size_t seed=0) const
      {
	for (/**/; first != last; ++ first)
	  seed = hasher_type::operator()(first->id(), seed);
	return seed;
      }
      
      template <typename Iterator, typename __Phrase>
      inline
      bool equal_phrase(Iterator first, Iterator last, const __Phrase& x) const
      {
	return x.size() == std::distance(first, last) && std::equal(first, last, x.begin());
      }

      
      double ngram_score(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge) const
      {
	const rule_type& rule = *(edge.rule);
	const rule_type::symbol_set_type& target = rule.target;
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	
	
	// we will reserve enough space so that buffer's memory will not be re-allocated.
	buffer.clear();
	buffer.reserve(target.size() + (order * 2) * states.size());
	
	double score = 0.0;
	
	int star_first = -1;
	int star_last = -1;
	
	buffer_type::const_iterator iter_first = buffer.begin();
	
	rule_type::symbol_set_type::const_iterator titer_end = target.end();
	for (rule_type::symbol_set_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer) {
	  if (titer->is_non_terminal()) {
	    const int antecedent_index = titer->non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      throw std::runtime_error("this is a non-terminal, but no index!");
	    
	    const symbol_type* context = reinterpret_cast<const symbol_type*>(states[antecedent_index]);
	    
	    buffer_type::const_iterator citer = buffer.end();
	    
	    buffer.insert(buffer.end(), context, std::find(context, context + order * 2, vocab_type::EMPTY));
	    
	    buffer_type::const_iterator citer_end = buffer.end();
	    for (/**/; citer != citer_end && *citer != vocab_type::STAR; ++ citer)
	      if ((citer + 1) - iter_first == order) {
		score += ngram.logprob(iter_first, citer + 1);
		++ iter_first;
	      }
	    
	    if (citer != citer_end) {
	      // update to the next of STAR
	      iter_first = citer + 1;
	      
	      if (star_first < 0)
		star_first = citer - buffer.begin();
	      star_last = citer - buffer.begin();
	    }
	    
	  } else if (*titer != vocab_type::EPSILON) {
	    // we will ignore epsilon...
	    buffer.push_back(*titer);
	    
	    const bool score_ngram = (buffer.end() - iter_first == order);
	    
	    if (score_ngram) {
	      buffer_type::const_iterator iter_end = buffer.end();
	      score += ngram.logprob(iter_first, iter_end);
	      ++ iter_first;
	    }
	  }
	}
	
	// construct state vector..
	symbol_type* context = reinterpret_cast<symbol_type*>(state);
	
	if (star_first >= 0) {
	  const int prefix_size = utils::bithack::min(star_first, order - 1);
	  const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), order - 1);
	  
	  std::copy(buffer.begin(), buffer.begin() + prefix_size, context);
	  context[prefix_size] = vocab_type::STAR;
	  std::copy(buffer.end() - suffix_size, buffer.end(), context + prefix_size + 1);
	  context[prefix_size + 1 + suffix_size] = vocab_type::EMPTY;
	} else {
	  if (buffer.size() <= order - 1) {
	    std::copy(buffer.begin(), buffer.end(), context);
	    context[buffer.size()] = vocab_type::EMPTY;
	  } else {
	    std::copy(buffer.begin(), buffer.begin() + order - 1, context);
	    context[order] = vocab_type::STAR;
	    std::copy(buffer.end() - (order - 1), buffer.end(), context + order);
	    context[order * 2 - 1] = vocab_type::EMPTY;
	  }
	}

	
	return score;
      }

      
      double ngram_final_score(state_ptr_type& state) const
      {
	const symbol_type* context      = reinterpret_cast<const symbol_type*>(state);
	const symbol_type* context_end  = std::find(context, context + order * 2, vocab_type::EMPTY);
	const symbol_type* context_star = std::find(context, context_end, vocab_type::STAR);
	
	buffer_id_type& buffer_id = const_cast<buffer_id_type&>(buffer_id_impl);
	
	buffer_id.clear();
	buffer_id.reserve(order * 2 + 2);
	
	buffer_id.push_back(ngram.index.vocab()[vocab_type::BOS]);
	
	double score = 0.0;
	for (const symbol_type* citer = context; citer != context_star; ++ citer) {
	  buffer_id.push_back(ngram.index.vocab()[*citer]);
	  score += ngram.logprob(std::max(buffer_id.begin(), buffer_id.end() - order), buffer_id.end());
	}

	if (context_star != context_end) {
	  buffer_id.clear();
	  for (const symbol_type* citer = context_star; citer != context_end; ++ citer)
	    buffer_id.push_back(ngram.index.vocab()[*citer]);
	}
	
	buffer_id.push_back(ngram.index.vocab()[vocab_type::EOS]);
	score += ngram.logprob(std::max(buffer_id.begin(), buffer_id.end() - order), buffer_id.end());
	
	return score;
	
#if 0
	const symbol_type* context = reinterpret_cast<const symbol_type*>(state);

	buffer_type& buffer       = const_cast<buffer_type&>(buffer_impl);
	
	buffer.clear();
	buffer.push_back(vocab_type::BOS);
	buffer.insert(buffer.end(), context, std::find(context, context + order * 2, vocab_type::EMPTY));
	buffer.push_back(vocab_type::EOS);

	buffer_type::const_iterator iter_begin = buffer.begin();
	buffer_type::const_iterator iter_end = buffer.end();
	
	buffer_type::const_iterator iter_first = iter_begin;
	
	double score = 0.0;
	
	// compute score until we reach STAR
	
	buffer_type::const_iterator iter = iter_begin + 1;
	for (/**/; iter != iter_end && *iter != vocab_type::STAR; ++ iter) {
	  score += ngram.logprob(iter_first, iter + 1);
	  iter_first += (iter + 1 - iter_first == order);
	}
	
	// we will rescore the EOS...
	if (iter != iter_end && iter + 1 != iter_end)
	  score += ngram.logprob(std::max(iter + 1, iter_end - order), iter_end);
	
	return score;
#endif
      }
      
      double ngram_estimate(state_ptr_type& state) const
      {
	const symbol_type* context = reinterpret_cast<const symbol_type*>(state);
	const symbol_type* context_end = context + order - 1;
	
	if (*context == vocab_type::EMPTY)
	  return 0.0;
	
	buffer_id_type& buffer_id = const_cast<buffer_id_type&>(buffer_id_impl);
	buffer_id.clear();
	
	double estimate = 0.0;
	for (const symbol_type* citer = context; citer != context_end && *citer != vocab_type::STAR && *citer != vocab_type::EMPTY; ++ citer) {
	  buffer_id.push_back(ngram.index.vocab()[*citer]);
	  bool estimated = false;
	  double logbound = ngram.logbound(buffer_id.begin(), buffer_id.end(), estimated);
	  
	  if (estimated && logbound < 0.0)
	    logbound *= decays[buffer_id.size()];
	  
	  estimate += logbound;
	}
	return estimate;
      }

      decay_set_type decays;
      
      // actual buffers...
      buffer_type    buffer_impl;
      buffer_id_type buffer_id_impl;
      
      // ngrams
      ngram_type     ngram;
      int            order;
    };
    
    
    NGram::NGram(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "ngram")
	throw std::runtime_error("is this really ngram feature function? " + parameter);

      parameter_type::const_iterator fiter = param.find("file");
      if (fiter == param.end() || ! boost::filesystem::exists(fiter->second))
	throw std::runtime_error("no ngram file?");

      const path_type path = fiter->second;
      
      int order = 3;
      parameter_type::const_iterator oiter = param.find("order");
      if (oiter != param.end())
	order = boost::lexical_cast<int>(oiter->second);
      if (order <= 0)
	throw std::runtime_error("invalid ngram order");
      
      std::auto_ptr<impl_type> ngram_impl(new impl_type(path, order));
      
      // two contexts (order - 1) for each edge, with two separator..
      base_type::__state_size = sizeof(symbol_type) * ngram_impl->order * 2;
      
      parameter_type::const_iterator niter = param.find("name");
      if (niter != param.end())
	base_type::__feature_name = niter->second;
      else
	base_type::__feature_name = std::string("ngram");
      
      pimpl = ngram_impl.release();
    }
    
    NGram::~NGram() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    void NGram::operator()(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge,
			   feature_set_type& features,
			   feature_set_type& estimates) const
    {
      features[base_type::feature_name()] = pimpl->ngram_score(state, states, edge);
      estimates[base_type::feature_name()] = pimpl->ngram_estimate(state);
    }
    
    void NGram::operator()(state_ptr_type& state,
			   feature_set_type& features) const
    {
      features[base_type::feature_name()] = pimpl->ngram_final_score(state);
    }
  };
};
