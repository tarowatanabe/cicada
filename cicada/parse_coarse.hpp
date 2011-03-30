// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_COARSE__HPP__
#define __CICADA__PARSE_COARSE__HPP__ 1

#include <vector>
#include <algorithm>
#include <set>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/b_heap.hpp>
#include <utils/std_heap.hpp>
#include <utils/bithack.hpp>
#include <utils/array_power2.hpp>


#include <google/dense_hash_map>
#include <google/dense_hash_set>

namespace cicada
{
  // coarse-to-fine parsing
  // input is a set of grammars, or use iterators
  // 
  // vector<grammar_type> grammar_set_type;
  //
  
  template <typename Semiring, typename Function>
  struct ParseCoarse
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;

    typedef Lattice    lattice_type;
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef attribute_set_type::attribute_type attribute_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    typedef Semiring semiring_type;
    typedef Semiring score_type;
    
    typedef Function function_type;
    
    typedef std::vector<grammar_type, std::allocator<grammar_type> > grammar_set_type;
    typedef std::vector<double, std::allocator<double> > threshold_set_type;

    class LabelScoreSet : public google::dense_hash_map<symbol_type, score_type, boost::hash<symbol_type>, std::equal_to<symbol_type> >
    {
    public:
      typedef google::dense_hash_map<symbol_type, score_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > label_score_set_type;
      
      LabelScoreSet() : label_score_set_type() { label_score_set_type::set_empty_key(symbol_type()); }
    };
    typedef LabelScoreSet label_score_set_type;
    typedef utils::chart<label_score_set_type, std::allocator<label_score_set_type > > label_score_chart_type;

    struct CoarseSymbol
    {
      CoarseSymbol(const int __bits) : bits(__bits) {}
      
      symbol_type operator()(const symbol_type& symbol)
      {
	if (! symbol.is_non_terminal()) return symbol;
	
	const size_t cache_pos = hash_value(symbol) & (caches.size() - 1);
	cache_type& cache = caches[cache_pos];
	if (cache.symbol != symbol) {
	  cache.symbol = symbol;
	  cache.coarse = symbol.coarse(bits);
	}
	return cache.coarse;
      }
      
      struct Cache
      {
	symbol_type symbol;
	symbol_type coarse;
	
	Cache() : symbol(), coarse() {}
      };
      typedef Cache cache_type;
      typedef utils::array_power2<cache_type, 1024 * 4, std::allocator<cache_type> > cache_set_type;
      
      cache_set_type caches;
      int bits;
    };
    
    struct CoarseSimple
    {
      symbol_type operator()(const symbol_type& symbol)
      {
	if (! symbol.is_non_terminal()) return symbol;
	
	const utils::piece piece = symbol.non_terminal_strip();
	
	// default X
	return (piece.find('^') != utils::piece::npos() ? "[x^]" : "[x]");
      }
    };
    
    template <typename Coarser>
    struct PruneCoarse
    {
      PruceCoarse(Coarser __coarser) : coarser(__coarser) {}
      
      bool operator()(const int first, const int last, const symbol_type& label)
      {
	const label_score_set_type& labels = prunes(first, last);
	
	label_score_set_type::const_iterator liter = labels.find(label);
	if (liter != labels.end())
	  return liter->second < cutoff;
	
	label_score_set_type::const_iterator citer = labels.find(coarser(label));
	return (citr == labels.end() || citer->second < cutoff);
      }
      
      const label_score_chart_type& prunes;
      const score_type cutoff;
      Coarser coarser;
    };
    
    struct PruneNone
    {
      bool operator()(const int first, const int last, const symbol_type& label) const { return false; }
    };
    
    struct ParseCKY
    {
      //
      // CKY parser... but we will not a construct hypergraph, but a tabular structure...
      //
      
      //
      // output is chart structure with label-score
      //
      
      ParseCKY(const symbol_type& __goal,
	       const grammar_type& __grammar,
	       const bool __yield_source=false,
	       const bool __treebank=false,
	       const bool __pos_mode=false)
	: goal(__goal), grammar(__grammar), threshold(__threshold), yield_source(__yield_source), treebank(__treebank), pos_mode(__pos_mode) {}
      
      template <typename Pruner>
      void operator()(const lattice_type& lattice, label_score_chart_type& scores, const Pruner& pruner)
      {
	
	
      }
      
      const symbol_type goal;
      const grammar_type& grammar;
      const bool yield_source;
      const bool treebank;
      const bool pos_mode;
    };
    
    template <typename IteratorGrammar, typename IteratorThreshold>
    ParseCoarse(const symbol_type& __goal,
		IteratorGrammar gfirst, IteratorGrammar glast,
		IteratorThreshold tfirst, IteratorThreshold tlast,
		const function_type& __function,
		const bool __yield_source=false,
		const bool __treebank=false,
		const bool __pos_mode=false)
      : goal(__goal),
	grammars(gfirst, glast),
	thresholds(tfirst, tlast),
	function(__function),
	yield_source(__yield_source),
	treebank(__treebank),
	pos_mode(__pos_mode),
	attr_span_first("span-first"),
	attr_span_last("span-last")
    {
      if (grammars.empty())
	throw std::runtime_error("no grammar?");
      if (thresholds.size() + 1 != grammar.size())
	throw std::runtime_error("do we have enough threshold parameters for grammars?");
    }
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      if (lattice.empty()) return;
      
      
      
    }
    
  private:
    const symbol_type goal;
    grammar_set_type   grammars;
    threshold_set_type thresholds;
    
    const function_type& function;

    const bool yield_source;
    const bool treebank;
    const bool pos_mode;
    const attribute_type attr_span_first;
    const attribute_type attr_span_last;
  };
  
  
  
  
};

#endif
