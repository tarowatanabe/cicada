#include <stdexcept>
#include <memory>

#include "cicada/ngram.hpp"
#include "cicada/feature_ngram.hpp"
#include "cicada/parameter.hpp"

namespace cicada
{
  namespace feature
  {
    class NGramImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Vocab  vocab_type;
      
      typedef cicada::NGram ngram_type;
      
      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;
      
      typedef feature_function_type::rule_type rule_type;
      
    public:
      typedef boost::filesystem::path path_type;
      
      NGramImpl(const path_type& __path, const int __order)
	: ngram(__path), order(__order)
      {
	order = utils::bithack::min(order, ngram.index.order());
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
	const symbol_type* context = reinterpret_cast<const symbol_type*>(state);

	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	
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
      }
      
      double ngram_estimate(state_ptr_type& state) const
      {
	const symbol_type* context = reinterpret_cast<const symbol_type*>(state);
	const symbol_type* context_end = context + order - 1;
	
	double score = 0.0;
	if (*context != vocab_type::EMPTY)
	  for (const symbol_type* citer = context + 1; citer != context_end && *citer != vocab_type::STAR && *citer != vocab_type::EMPTY; ++ citer)
	    score += ngram.logbound(context, citer + 1);
	
	return score;
      }
      
      ngram_type  ngram;
      buffer_type buffer_impl;
      int         order;
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
