// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MBR__HPP__
#define __CICADA__MBR__HPP__ 1

#include <cicada/kbest.hpp>
#include <cicada/expected_ngram.hpp>

#include <utils/sgi_hash_map.hpp>

namespace cicada
{
  
  template <typename Traversal, typename Function, typename Filter, typename Extract>
  struct MBR
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::id_type   id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef rule_type::vocab_type      vocab_type;
    typedef rule_type::symbol_type     symbol_type;
    typedef rule_type::symbol_set_type symbol_set_type;
    
    typedef Traversal traversal_type;
    typedef Function  function_type;
    typedef Filter    filter_type;
    typedef Extract   extract_type;
    
    typedef typename traversal_type::value_type yield_type;
    typedef typename function_type::value_type  semiring_type;
    typedef typename function_type::value_type  weight_type;
    typedef typename extract_type::value_type   sentence_type;
    
    typedef symbol_type     word_type;
    typedef symbol_set_type ngram_type;

    
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<ngram_type, double, boost::hash<ngram_type>, std::equal_to<ngram_type>, std::allocator<std::pair<const ngram_type, double> > > ngram_set_type;
#else
    typedef sgi::hash_map<ngram_type, double, boost::hash<ngram_type>, std::equal_to<ngram_type>, std::allocator<std::pair<const ngram_type, double> > > ngram_set_type;
#endif

    typedef KBest<Traversal, Function, Filter> kbest_derivation_type;
    typedef typename kbest_type::derivation_derivation_type derivation_type;

    typedef std::pair<weight_type, yield_type> weight_yield_type;
    typedef std::vector<weight_yield_type, std::allocator<weight_yield_type> > weight_yield_set_type;
    
    
    MBR(const hypergraph_type& __graph,
	const int __order,
	const int __kbest,
	const traversal_type& __traversal,
	const function_type& __function,
	const filter_type& __filter,
	const extract_type& __extract,
	const bool __yield_source=false)
      : graph(__graph),
	order(__order),
	kbest(__kbest)
	yield_source(__yield_source)
    {
      // compute expected ngram counts
      ngram_set_type ngrams;
      
      cicada::expected_ngram(graph, __function, ngrams, order, false, yield_source);
      
      // compute kbest-translations
      kbest_derivation_type derivations(__graph, __kbest, __traversal, __function, __filter);

      weight_yields.clear();
      weight_yields.resize(kbest);
      
      for (int k = 0; k < kbest; ++ k)
	if (! derivations(k, weight_yields[k].second, weight_yields[k].first)) {
	  weight_yields.resize(k);
	  break;
	}
      
      // now, we will perform additional compuation..
      double expected_length = 0.0;
      {
	ngram_set_type::const_iterator niter_end = ngrams.end();
	for (ngram_set_type::const_iterator niter = ngrams.begin(); niter != niter_end; ++ niter)
	  if (niter->first.size() == 1)
	    expected_length += niter->second;
      }
      
      weight_yield_set_type::iterator yiter_end = weight_yields.end();
      for (weight_yield_set_type::iterator yiter = weight_yields.begin(); yiter != yiter_end; ++ yiter) {
	const sentence_type& sentence = __extract(yiter->second);
	
	const double pnalty = std::min(0.0, 1.0 - (expected_length / sentence.size()));
	const double factor = 1.0 / order;
	
	double score = penalty;
	for (int n = 1; n <= order; ++ n) {
	  ngram_set_type ngrams_local;
	  
	  sentence_type::const_iterator siter_end = sentence.end();
	  for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	    // ngram at [iter, iter + n)
	    for (sentence_type::const_iterator iter = siter; iter + n < siter_end; ++ iter)
	      ngrams_local[ngram_type(iter, iter + 2)] += 1;
	  }
	  
	  double count = 0.0;
	  ngram_set_type::const_iterator niter_end = ngrams_local.end();
	  for (ngram_set_type::const_iterator niter = ngrams_local.begin(); niter != niter_end; ++ niter) {
	    ngram_set_type::const_iterator niter_prime = ngrams.find(niter->first);
	    if (niter_prime != ngrams.end())
	      count += std::min(niter->second, niter_prime->second);
	  }
	  
	  if (count != 0.0)
	    score += factor * std::log(count / (sentence.size() - n + 1));
	}
	
	yiter->first = cicada::semiring::traits<weight_type>::log(std::exp(score));
      }
      
      std::sort(weight_yields.begin(), weight_yields.end());
    }
    
    bool operator()(int k, yield_type& yield, weight_type& weight)
    {
      if (k < weight_yields.size()) {
	yield  = weight_yields[k].second;
	weight = weight_yields[k].first;
	return true;
      } else {
	yield  = yield_type();
	weight = weight_type();
	return false;
      }
    }
    
  private:
    const hypergraph_type& graph;
    
    weight_yield_set_type weight_yields;
    
    const int  order;
    const int  kbest;
    const bool yield_source;
  };
  
};

#endif
