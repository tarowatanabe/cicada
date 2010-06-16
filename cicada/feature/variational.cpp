#include <stdexcept>
#include <memory>

#include "cicada/feature/variational.hpp"
#include "cicada/parameter.hpp"

#include "utils/compact_trie.hpp"

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

      typedef double logprob_type;
      
      typedef utils::compact_trie<symbol_type, logprob_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				  std::allocator<std::pair<const symbol_type, logprob_type> > > ngram_set_type;

      typedef std::vector<feature_type, std::allocator<feature_type> > feature_name_set_type;

      typedef cicada::FeatureFunction feature_function_type;

      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;

      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;

      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;

      
    public:
      VariationalImpl(const int __order)
	: order(__order)
      {
	feature_names.resize(order);
	for (int n = 0; n < order; ++ n)
	  feature_names[n] = "variational" + boost::lexical_cast<std::string>(n + 1);
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
	const phrase_type& target = rule.target;

	phrase_type::const_iterator titer_begin = target.begin();
	phrase_type::const_iterator titer_end   = target.end();
	
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
	  
	  ngram_score(buffer.begin(), buffer.end(), features);
	  
	  if (buffer.size() <= context_size)
	    std::copy(buffer.begin(), buffer.end(), context);
	  else {
	    std::copy(buffer.begin(), buffer.begin() + context_size, context);
	    context[context_size] = vocab_type::STAR;
	    std::copy(buffer.end() - context_size, buffer.end(), context + order);
	  }
	} else {
	  phrase_span_set_type& phrase_spans = const_cast<phrase_span_set_type&>(phrase_spans_impl);
	  phrase_spans.clear();
	  target.terminals(std::back_inserter(phrase_spans));
	
	  int star_first = -1;
	  int star_last  = -1;
	  
	  for (phrase_type::const_iterator titer = phrase_spans.front().first; titer != phrase_spans.front().second; ++ titer)
	    if (*titer != vocab_type::EPSILON)
	      buffer.push_back(*titer);
	  
	  ngram_score(buffer.begin(), buffer.end(), features);
	  
	  buffer_type::const_iterator biter_first = buffer.begin();
	  
	  phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	  for (phrase_span_set_type::const_iterator siter = phrase_spans.begin() + 1; siter != siter_end; ++ siter) {
	    const phrase_span_type& span = *siter;
	    
	    const int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      throw std::runtime_error("this is a non-terminal, but no index!");
	    
	    const symbol_type* context = reinterpret_cast<const symbol_type*>(states[antecedent_index]);
	    const symbol_type* context_end  = std::find(context, context + order * 2, vocab_type::EMPTY);
	    const symbol_type* context_star = std::find(context, context_end, vocab_type::STAR);
	    
	    buffer_type::const_iterator biter = buffer.end();
	    
	    buffer.insert(buffer.end(), context, context_star);
	    
	    buffer_type::const_iterator biter_end = buffer.end();

	    ngram_score(biter_first, biter, biter_end, features);
	    
	    // insert star!
	    if (context_star != context_end) {
	      
	      biter_first = buffer.end() + 1;
	      
	      star_last = buffer.size() + 1;
	      if (star_first < 0)
		star_first = buffer.size() + 1;
	      
	      buffer.insert(buffer.end(), context_star, context_end);
	    }
	    
	    // use of this edge's terminals
	    {
	      buffer_type::const_iterator biter = buffer.end();
	      
	      buffer.insert(buffer.end(), span.first, span.second);
	      
	      buffer_type::const_iterator biter_end = buffer.end();
	      
	      ngram_score(biter_first, biter, biter_end, features);
	    }
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
	  } else {
	    if (buffer.size() <= order - 1)
	      std::copy(buffer.begin(), buffer.end(), context);
	    else {
	      std::copy(buffer.begin(), buffer.begin() + context_size, context);
	      context[context_size] = vocab_type::STAR;
	      std::copy(buffer.end() - context_size, buffer.end(), context + order);
	    }
	  }
	}
      }

      void ngram_score(const state_ptr_type& state,
		       feature_set_type& features) const
      {
	const symbol_type* context      = reinterpret_cast<const symbol_type*>(state);
	const symbol_type* context_end  = std::find(context, context + order * 2, vocab_type::EMPTY);
	const symbol_type* context_star = std::find(context, context_end, vocab_type::STAR);
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	buffer.clear();
	
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
      phrase_span_set_type phrase_spans_impl;

      feature_name_set_type feature_names;
      int order;
    };
    
    Variational::Variational(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (param.name() != "variational")
	throw std::runtime_error("is this variational feature?" + parameter);
      
      int order = 3;
      parameter_type::const_iterator oiter = param.find("order");
      if (oiter != param.end())
	order = boost::lexical_cast<int>(oiter->second);
      if (order <= 0)
	throw std::runtime_error("invalid ngram order");
      
      base_type::__state_size = sizeof(symbol_type) * order * 2;
      
      pimpl = new impl_type(order);
    }
    
    Variational::~Variational() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    void Variational::operator()(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates) const
    {
      pimpl->ngram_score(state, states, edge, features);
    }
    
    void Variational::operator()(const state_ptr_type& state,
				 feature_set_type& features) const
    {
      pimpl->ngram_score(state, features);
    }
    
    void Variational::clear()
    {
      pimpl->clear();
    }

      
    void Variational::insert(const ngram_type& ngram, const double& logprob)
    {
      pimpl->insert(ngram.begin(), ngram.end(), logprob);
    }
    
  };
};
