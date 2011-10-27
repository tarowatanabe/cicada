//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <memory>

#include "cicada/feature/variational.hpp"
#include "cicada/parameter.hpp"
#include "cicada/inside_outside.hpp"
#include "cicada/semiring.hpp"

#include "utils/compact_trie_dense.hpp"
#include "utils/sgi_hash_map.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"

#include <boost/lexical_cast.hpp>

namespace cicada
{
  namespace feature
  {
    class VariationalImpl
    {
      public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Symbol   word_type;
      typedef cicada::Feature  feature_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;

      typedef HyperGraph hypergraph_type;
      
      typedef double logprob_type;
      
      typedef utils::compact_trie_dense<symbol_type, logprob_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
					std::allocator<std::pair<const symbol_type, logprob_type> > > ngram_set_type;

      typedef std::vector<feature_type, std::allocator<feature_type> > feature_name_set_type;

      typedef cicada::FeatureFunction feature_function_type;

      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;

      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;

      
    public:
      VariationalImpl(const int __order)
	: ngrams(symbol_type()), order(__order), no_bos_eos(false)
      {
	feature_names.resize(order);
	for (int n = 0; n < order; ++ n)
	  feature_names[n] = "variational" + utils::lexical_cast<std::string>(n + 1);
      }
      
    public:

      template <typename Iterator>
      void ngram_score(Iterator first, Iterator iter, Iterator last, feature_set_type& features) const
      {
	const int context_size = order - 1;
	
	first = std::max(iter - context_size, first);
	
	for (/**/; first != iter; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter2 = first; iter2 != std::min(first + order, last); ++ iter2) {
	    id = ngrams.find(id, *iter2);
	    
	    if (ngrams.is_root(id)) break;
	    if (iter2 < iter) continue;
	    
	    features[feature_names[iter2 - first]] += ngrams[id];
	  }
	}
      }

      template <typename Iterator>
      void ngram_score(Iterator first, Iterator last, feature_set_type& features) const
      {
	for (/**/; first != last; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter = first; iter != std::min(first + order, last); ++ iter) {
	    id = ngrams.find(id, *iter);
	    
	    if (ngrams.is_root(id)) break;
	    
	    features[feature_names[iter - first]] += ngrams[id];
	  }
	}
      }

      void ngram_score(state_ptr_type& state,
		       const state_ptr_set_type& states,
		       const edge_type& edge,
		       feature_set_type& features) const
      {
	const int context_size = order - 1;
	const rule_type& rule = *(edge.rule);
	const phrase_type& phrase = rule.rhs;

	phrase_type::const_iterator titer_begin = phrase.begin();
	phrase_type::const_iterator titer_end   = phrase.end();
	
	// we will reserve enough space so that buffer's memory will not be re-allocated.
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	buffer.reserve((titer_end - titer_begin) + (order * 2) * states.size());
	
	if (states.empty()) {
	  symbol_type* context = reinterpret_cast<symbol_type*>(state);
	  std::fill(context, context + order * 2, vocab_type::EMPTY);
	  
	  for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	    if (*titer != vocab_type::EPSILON)
	      buffer.push_back(*titer);
	  
	  if (! ngrams.empty())
	    ngram_score(buffer.begin(), buffer.end(), features);
	  
	  if (static_cast<int>(buffer.size()) <= context_size)
	    std::copy(buffer.begin(), buffer.end(), context);
	  else {
	    std::copy(buffer.begin(), buffer.begin() + context_size, context);
	    context[context_size] = vocab_type::STAR;
	    std::copy(buffer.end() - context_size, buffer.end(), context + order);
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
	      
	      const symbol_type* context = reinterpret_cast<const symbol_type*>(states[antecedent_index]);
	      const symbol_type* context_end  = std::find(context, context + order * 2, vocab_type::EMPTY);
	      const symbol_type* context_star = std::find(context, context_end, vocab_type::STAR);
	      
	      if (biter != buffer.end()) {
		if (! ngrams.empty()) {
		  if (biter_first != biter)
		    ngram_score(biter_first, biter, buffer.end(), features);
		  ngram_score(biter, buffer.end(), features);
		}
		biter = buffer.end();
	      }
	      
	      buffer.insert(buffer.end(), context, context_star);
	      if (! ngrams.empty())
		if (biter_first != biter && biter != buffer.end())
		  ngram_score(biter_first, biter, buffer.end(), features);
	      biter = buffer.end();
	      
	      if (context_star != context_end) {
		star_first = utils::bithack::branch(star_first < 0, static_cast<int>(buffer.size()) + 1, star_first);
		star_last  = buffer.size() + 1;
		
		biter_first = buffer.end() + 1;
		buffer.insert(buffer.end(), context_star, context_end);
		biter = buffer.end();
	      }
	    } else if (*titer != vocab_type::EPSILON)
	      buffer.push_back(*titer);
	  }
	  
	  if (biter != buffer.end()) {
	    if (! ngrams.empty()) {
	      if (biter_first != biter)
		ngram_score(biter_first, biter, buffer.end(), features);
	      ngram_score(biter, buffer.end(), features);
	    }
	    biter = buffer.end();
	  }
	  
	  // construct state vector..
	  symbol_type* context = reinterpret_cast<symbol_type*>(state);
	  std::fill(context, context + order * 2, vocab_type::EMPTY);
	  
	  if (star_first >= 0) {
	    const int prefix_size = utils::bithack::min(star_first, context_size);
	    const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
	    
	    std::copy(buffer.begin(), buffer.begin() + prefix_size, context);
	    context[prefix_size] = vocab_type::STAR;
	    std::copy(buffer.end() - suffix_size, buffer.end(), context + prefix_size + 1);
	  } else if (static_cast<int>(buffer.size()) <= context_size)
	    std::copy(buffer.begin(), buffer.end(), context);
	  else {
	    std::copy(buffer.begin(), buffer.begin() + context_size, context);
	    context[context_size] = vocab_type::STAR;
	    std::copy(buffer.end() - context_size, buffer.end(), context + order);
	  }
	}
      }

      void ngram_score(const state_ptr_type& state,
		       feature_set_type& features) const
      {
	if (ngrams.empty())
	  return;

	const symbol_type* context      = reinterpret_cast<const symbol_type*>(state);
	const symbol_type* context_end  = std::find(context, context + order * 2, vocab_type::EMPTY);
	const symbol_type* context_star = std::find(context, context_end, vocab_type::STAR);
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	
	if (no_bos_eos) {
	  buffer.insert(buffer.end(), context, context_star);
	  
	  if (! buffer.empty())
	    ngram_score(buffer.begin(), buffer.begin() + (buffer.front() == vocab_type::BOS), buffer.end(), features);
	} else {
	  buffer.push_back(vocab_type::BOS);
	  buffer.insert(buffer.end(), context, context_star);
	  
	  ngram_score(buffer.begin(), buffer.begin() + 1, buffer.end(), features);
	  
	  if (context_star != context_end) {
	    buffer.clear();
	    buffer.insert(buffer.end(), context_star + 1, context_end);
	  }
	  buffer.push_back(vocab_type::EOS);
	  
	  ngram_score(buffer.begin(), buffer.end() - 1, buffer.end(), features);
	}
      }
      
      void clear()
      {
	ngrams.clear();
      }
      
      template <typename Iterator>
      void insert(Iterator first, Iterator last, const logprob_type& logprob)
      {
	ngrams[ngrams.insert(first, last)] =  logprob;
      }
      
    private:
      ngram_set_type ngrams;

      buffer_type          buffer_impl;

      feature_name_set_type feature_names;
      
    public:
      int order;
      bool no_bos_eos;
    };
    
    Variational::Variational(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "variational")
	throw std::runtime_error("is this variational feature?" + parameter);
      
      int order = 3;
      bool no_bos_eos = false;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "no-bos-eos")
	  no_bos_eos = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for variational: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (order <= 0)
	throw std::runtime_error("invalid ngram order");
      
      base_type::__state_size = sizeof(symbol_type) * order * 2;
      base_type::__feature_name = "variational";
      
      pimpl = new impl_type(order);
      pimpl->no_bos_eos = no_bos_eos;
    }
    
    Variational::~Variational() { std::auto_ptr<impl_type> tmp(pimpl); }

    Variational::Variational(const Variational& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}
    
    Variational& Variational::operator=(const Variational& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }

    void Variational::apply(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features,
			    const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));

      pimpl->ngram_score(state, states, edge, features);

      if (final)
	pimpl->ngram_score(state, features);
    }

    void Variational::apply_coarse(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   const bool final) const
    {
    }
    void Variational::apply_predict(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    const bool final) const
    {}
    void Variational::apply_scan(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 const int dot,
				 feature_set_type& features,
				 const bool final) const
    {}
    void Variational::apply_complete(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     feature_set_type& features,
				     const bool final) const
    {
      apply(state, states, edge, features, final);
    }

    
    void Variational::clear()
    {
      pimpl->clear();
    }

    void Variational::assign(const size_type& id,
			     const hypergraph_type& hypergraph,
			     const lattice_type& lattice,
			     const span_set_type& spans,
			     const sentence_set_type& targets,
			     const ngram_count_set_type& ngram_counts)
    {
      typedef ngram_count_set_type::ngram_type ngram_type;

      pimpl->clear();
      
      ngram_count_set_type history;
      
      ngram_count_set_type::const_iterator niter_end = ngram_counts.end();
      for (ngram_count_set_type::const_iterator niter = ngram_counts.begin(); niter != niter_end; ++ niter)
	history[ngram_type(niter->first.begin(), niter->first.end() - 1)] += niter->second;
      
      for (ngram_count_set_type::const_iterator niter = ngram_counts.begin(); niter != niter_end; ++ niter)
	pimpl->insert(niter->first.begin(), niter->first.end(),
		      std::log(niter->second / history[ngram_type(niter->first.begin(), niter->first.end() - 1)]));
    };
  };
};
