// -*- mode: c++ -*-

#ifndef __CICADA__MBR__HPP__
#define __CICADA__MBR__HPP__ 1

#include <cicada/kbest.hpp>
#include <cicada/expected_ngram.hpp>

#include <utils/sgi_hash_map.hpp>

namespace cicada
{
  
  template <typename Traversal, typename Function, typename Filter>
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
    
    typedef typename traversal_type::value_type yield_type;
    typedef typename function_type::value_type  semiring_type;
    typedef typename function_type::value_type  weight_type;
    
    typedef symbol_type     word_type;
    typedef symbol_set_type ngram_type;

    
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<ngram_type, double, boost::hash<ngram_type>, std::equal_to<ngram_type>, std::allocator<std::pair<const ngram_type, double> > > ngram_set_type;
#else
    typedef sgi::hash_map<ngram_type, double, boost::hash<ngram_type>, std::equal_to<ngram_type>, std::allocator<std::pair<const ngram_type, double> > > ngram_set_type;
#endif

    typedef KBest<Traversal, Function, Filter> kbest_derivation_type;
    typedef typename kbest_type::derivation_derivation_type derivation_type;
    
    
    MBR(const hypergraph_type& __graph,
	const int __order,
	const int __kbest,
	const traversal_type& __traversal,
	const function_type& __function,
	const filter_type& __filter,
	const bool __yield_source=false)
      : graph(__graph),
	kbest_derivation(graph, __kbest, __traversal, __function, __filter),
	order(__order),
	kbest(__kbest)
	yield_source(__yield_source) {}
    
    bool operator()(int k, yield_type& yield, weight_type& weight)
    {
      
      
      
    }
    
  private:
    
    template <typename YieldSet, typename ExtractYield, typename ExtractScore>
    void operator()(const hypergraph_type& graph, YieldSet& yields, ExtractYield extract_yield, ExtractScore extract_score)
    {
      ngram_set_type ngrams;
      
      cicada::expected_ngram(garph, function, ngrams, order, false, yield_source);
      
      
      
    }

    const hypergraph_type& graph;

    kbest_derivation_type kbest_derivation;
    
    const function_type function;
    
    const int  order;
    const int  kbest;
    const bool yield_source;
  };
  
};

#endif
