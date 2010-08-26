#include <stdexcept>
#include <memory>

#include "cicada/feature/variational.hpp"
#include "cicada/parameter.hpp"
#include "cicada/inside_outside.hpp"
#include "cicada/semiring.hpp"

#include "utils/compact_trie_dense.hpp"
#include "utils/sgi_hash_map.hpp"

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

      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;

      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;

      
    public:
      VariationalImpl(const int __order)
	: ngrams(symbol_type()), order(__order)
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
	  
	  if (! ngrams.empty())
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
	  
	  if (! ngrams.empty())
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
	    
	    if (! ngrams.empty())
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
	      
	      if (! ngrams.empty())
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
	    if (buffer.size() <= context_size)
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
	if (ngrams.empty())
	  return;

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
      
    public:
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
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "order") == 0)
	  order = boost::lexical_cast<int>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for variational: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (order <= 0)
	throw std::runtime_error("invalid ngram order");
      
      base_type::__state_size = sizeof(symbol_type) * order * 2;
      base_type::__feature_name = "variational";
      
      pimpl = new impl_type(order);
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
    
    template <typename FeaturePrefix, typename Feature>
    inline
    bool equal_prefix(const FeaturePrefix& prefix, const Feature& x)
    {
      return x.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), x.begin());
    }

    void Variational::operator()(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {
      const std::string& __feature_prefix = base_type::feature_name();
      for (feature_set_type::iterator fiter = features.begin(); fiter != features.end(); /**/)
	if (equal_prefix(__feature_prefix, fiter->first))
	  features.erase(fiter ++);
	else
	  ++ fiter;

      pimpl->ngram_score(state, states, edge, features);

      if (final)
	pimpl->ngram_score(state, features);
    }

    int Variational::order() const
    {
      return pimpl->order;
    }
    
    void Variational::clear()
    {
      pimpl->clear();
    }

    struct variational_function
    {
      typedef cicada::semiring::Logprob<double> weight_type;
      typedef weight_type value_type;

      typedef Variational::weight_set_type weight_set_type;

      const weight_set_type& weights;

      variational_function(const weight_set_type& __weights)
	: weights(__weights) {}

      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
      }
    };

    template <typename Iterator, typename Counts>
    void collect_counts(Iterator first, Iterator iter, Iterator last, const double& weight, Counts& counts, const int order)
    {
      const int context_size = order - 1;
      
      first = std::max(iter - context_size, first);
	
      for (/**/; first != iter; ++ first)
	for (Iterator iter2 = iter; iter2 != std::min(first + order, last); ++ iter2)
	  counts[typename Counts::key_type(first, iter2 + 1)] += weight;
    }

    template <typename Iterator, typename Counts>
    void collect_counts(Iterator first, Iterator last, const double& weight, Counts& counts, const int order)
    {
      for (/**/; first != last; ++ first)
	for (Iterator iter = first; iter != std::min(first + order, last); ++ iter)
	  counts[typename Counts::key_type(first, iter + 1)] += weight;
    }
    
      
    void Variational::insert(const hypergraph_type& graph, const weight_set_type& weights)
    {
      typedef variational_function::weight_type weight_type;
      
      typedef rule_type::symbol_set_type ngram_type;
      typedef rule_type::symbol_set_type phrase_type;

#ifdef HAVE_TR1_UNORDERED_MAP
      typedef std::tr1::unordered_map<ngram_type, double, boost::hash<ngram_type>, std::equal_to<ngram_type>,
	std::allocator<std::pair<ngram_type, double> > > ngram_set_type;
#else
      typedef sgi::hash_map<ngram_type, double, boost::hash<ngram_type>, std::equal_to<ngram_type>,
	std::allocator<std::pair<ngram_type, double> > > ngram_set_type;
#endif

      typedef std::pair<ngram_type, ngram_type> context_type;

      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;

      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;

      pimpl->clear();
      if (! graph.is_valid())
	return;
      
      const int order = pimpl->order;
      const int context_size = order - 1;

      ngram_set_type ngrams;

      phrase_span_set_type phrase_spans;
      buffer_type buffer;
      
      std::vector<weight_type, std::allocator<weight_type> > inside(graph.nodes.size());
      std::vector<weight_type, std::allocator<weight_type> > outside(graph.nodes.size());
      std::vector<context_type, std::allocator<context_type> > contexts(graph.nodes.size());
      
      cicada::inside(graph, inside, variational_function(weights));
      cicada::outside(graph, inside, outside, variational_function(weights));
      
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
	const node_type& node = *niter;

	const bool is_goal = (node.id == graph.goal);
	
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const edge_type& edge = graph.edges[*eiter];
	  
	  weight_type weight = cicada::semiring::traits<weight_type>::log(edge.features.dot(weights)) * inside[node.id] / inside[graph.goal];
	  edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	  for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	    weight *= outside[*niter];
	  
	  // we will construct ngram....
	  const phrase_type& target = edge.rule->target;

	  buffer.clear();
	  
	  if (edge.tails.empty()) {
	    phrase_type::const_iterator titer_end = target.end();
	    for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	      if (*titer != vocab_type::EPSILON)
		buffer.push_back(*titer);

	    collect_counts(buffer.begin(), buffer.end(), weight, ngrams, order);
	    
	    if (contexts[node.id].first.empty()) {
	      if (buffer.size() <= context_size)
		contexts[node.id] = std::make_pair(ngram_type(buffer.begin(), buffer.end()), ngram_type());
	      else
		contexts[node.id] = std::make_pair(ngram_type(buffer.begin(), buffer.begin() + context_size),
						   ngram_type(buffer.end() - context_size, buffer.end()));
	    }

	    if (is_goal) {
	      buffer.insert(buffer.begin(), vocab_type::BOS);
	      buffer.insert(buffer.end(), vocab_type::EOS);
	      
	      collect_counts(buffer.begin(), buffer.begin() + 1, buffer.end(), weight, ngrams, order);
	      collect_counts(buffer.begin(), buffer.end() - 1, buffer.end(), weight, ngrams, order);
	    }

	  } else {
	    buffer.reserve(target.size() + edge.tails.size() * order * 2);

	    phrase_spans.clear();
	    target.terminals(std::back_inserter(phrase_spans));	    
	    
	    int star_first = -1;
	    int star_last  = -1;
	    
	    for (phrase_type::const_iterator titer = phrase_spans.front().first; titer != phrase_spans.front().second; ++ titer)
	      if (*titer != vocab_type::EPSILON)
		buffer.push_back(*titer);
	    
	    collect_counts(buffer.begin(), buffer.end(), weight, ngrams, order);
	    
	    buffer_type::const_iterator biter_first = buffer.begin();
	    
	    phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	    for (phrase_span_set_type::const_iterator siter = phrase_spans.begin() + 1; siter != siter_end; ++ siter) {
	      const phrase_span_type& span = *siter;
	      
	      const int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	      if (antecedent_index < 0)
		throw std::runtime_error("this is a non-terminal, but no index!");

	      const context_type& context = contexts[edge.tails[antecedent_index]];
	      
	      buffer_type::const_iterator biter = buffer.end();
	      
	      buffer.insert(buffer.end(), context.first.begin(), context.first.end());
	      
	      buffer_type::const_iterator biter_end = buffer.end();
	      
	      collect_counts(biter_first, biter, biter_end, weight, ngrams, order);
	      
	      if (! context.second.empty()) {
		biter_first = buffer.end();
		
		star_last = buffer.size();
		if (star_first < 0)
		  star_first = buffer.size();
		
		buffer.insert(buffer.end(), context.second.begin(), context.second.end());
	      }
	      
	      {
		buffer_type::const_iterator biter = buffer.end();
	      
		buffer.insert(buffer.end(), span.first, span.second);
		
		buffer_type::const_iterator biter_end = buffer.end();
		
		collect_counts(biter_first, biter, biter_end, weight, ngrams, order);
	      }
	    }
	    
	    // setup history..
	    if (contexts[node.id].first.empty()) {
	      
	      if (star_first >= 0) {
		const int prefix_size = utils::bithack::min(star_first, context_size);
		const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
		
		contexts[node.id]= std::make_pair(ngram_type(buffer.begin(), buffer.begin() + prefix_size),
						  ngram_type(buffer.end() - suffix_size, buffer.end()));
		
	      } else {
		if (buffer.size() <= context_size)
		  contexts[node.id] = std::make_pair(ngram_type(buffer.begin(), buffer.end()), ngram_type());
		else
		  contexts[node.id] = std::make_pair(ngram_type(buffer.begin(), buffer.begin() + context_size),
						     ngram_type(buffer.end() - context_size, buffer.end()));
	      }
	    }
	    
	    if (is_goal) {
	      if (star_first >= 0) {
		const int prefix_size = utils::bithack::min(star_first, context_size);
		const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
		
		buffer.insert(buffer.begin(), vocab_type::BOS);
		buffer.insert(buffer.end(), vocab_type::EOS);

		collect_counts(buffer.begin(), buffer.begin() + 1, buffer.begin() + 1 + prefix_size, weight, ngrams, order);
		collect_counts(buffer.end() - suffix_size - 1, buffer.end() - 1, buffer.end(), weight, ngrams, order);
	      } else {
		buffer.insert(buffer.begin(), vocab_type::BOS);
		buffer.insert(buffer.end(), vocab_type::EOS);

		collect_counts(buffer.begin(), buffer.begin() + 1, buffer.end(), weight, ngrams, order);
		collect_counts(buffer.begin(), buffer.end() - 1, buffer.end(), weight, ngrams, order);
	      }
	    }
	  }
	}
      }

      ngram_set_type ngrams_history;
      {
	ngram_set_type::const_iterator niter_end = ngrams.end();
	for (ngram_set_type::const_iterator niter = ngrams.begin(); niter != niter_end; ++ niter)
	  ngrams_history[ngram_type(niter->first.begin(), niter->first.end() - 1)] += niter->second;
	
	for (ngram_set_type::const_iterator niter = ngrams.begin(); niter != niter_end; ++ niter)
	  pimpl->insert(niter->first.begin(), niter->first.end(), std::log(niter->second / ngrams_history[ngram_type(niter->first.begin(), niter->first.end() - 1)]));
      }
    }
    
  };
};
