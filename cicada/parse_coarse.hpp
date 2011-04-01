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
#include <cicada/span_vector.hpp>
#include <cicada/compose_cky.hpp>

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
  // Actually, we need different implementation....
  // How to handle this in the same framework.... separate this into different project...?
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
      
      symbol_type operator()(const symbol_type& symbol) const
      {
	if (! symbol.is_non_terminal()) return symbol;
	
	const size_t cache_pos = hash_value(symbol) & (caches.size() - 1);
	cache_type& cache = const_cast<cache_type&>(caches[cache_pos]);
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
      symbol_type operator()(const symbol_type& symbol) const
      {
	if (! symbol.is_non_terminal()) return symbol;
	
	return (symbol.binarized() ? "[x^]" : "[x]");
      }
    };
    
    template <typename Coarser>
    struct PruneCoarse
    {
      PruneCoarse(const label_score_chart_type& __prunes,
		  const score_type& __cutoff,
		  Coarser __coarser)
	: prunes(__prunes),
	  cutoff(__cutoff),
	  coarser(__coarser) {}

      bool operator()(const int first, const int last) const
      {
	return prunes(first, last).empty();
      }
      
      bool operator()(const int first, const int last, const symbol_type& label) const
      {
	if (prunes(first, last).empty()) return true;
	
	const label_score_set_type& labels = prunes(first, last);
	
	typename label_score_set_type::const_iterator liter = labels.find(label);
	if (liter != labels.end())
	  return liter->second < cutoff;
	
	typename label_score_set_type::const_iterator citer = labels.find(coarser(label));
	return (citer == labels.end() || citer->second < cutoff);
      }
      
      const label_score_chart_type& prunes;
      const score_type cutoff;
      Coarser coarser;
    };
    
    struct PruneNone
    {
      bool operator()(const int first, const int last) const
      {
	return false;
      }

      bool operator()(const int first, const int last, const symbol_type& label) const
      {
	return false;
      }
    };
    
    struct ParseCKY
    {
      //
      // CKY parser... but we will not a construct hypergraph, but a tabular structure...
      //
      // we keep two-passives, passives for keeping unary chain result (including zero-unary!)
      // and passives for keeping after non-unary...
      //
      
      typedef cicada::SpanVector::span_type span_type;
      typedef cicada::SpanVector::span_type tail_type;
      typedef utils::simple_vector<tail_type, std::allocator<tail_type> > tail_set_type;

      struct Edge
      {
	tail_set_type tails;
	score_type score;
	
	Edge() : tails(), score(cicada::semiring::traits<score_type>::one()) {}
	Edge(const tail_set_type& __tails, const score_type& __score) : tails(__tails), score(__score) {}
      };
      
      typedef Edge edge_type;
      typedef std::vector<edge_type, std::allocator<edge_type> > edge_set_type;
      
      struct Active
      {
	Active(const transducer_type::id_type& __node)
	  : node(__node), edge() {}
	Active(const transducer_type::id_type& __node,
	       const edge_type& __edge)
	  : node(__node), edge(__edge) {}
	
	transducer_type::id_type node;
	edge_type edge;
      };
      
      struct Passive
      {
	Passive(const span_type& __span) : span(__span), edges() {}
	
	span_type     span;
	edge_set_type edges;
      };
      
      struct Unary
      {
	symbol_type label;
	score_type  score;
	
	Unary() : label(), score(cicada::semiring::traits<score_type>::one()) {}
	Unary(const symbol_type& __label, const score_type& __score)
	  : label(__label), score(__score) {}
	template <typename Label, typename Score>
	Unary(const std::pair<Label, Score>& x)
	  : label(x.first), score(x.second) {}
      };
      
      typedef Active  active_type;
      typedef Passive passive_type;
      typedef Unary   unary_type;
      
      typedef utils::chunk_vector<active_type, 4096 / sizeof(active_type), std::allocator<active_type> > active_set_type;
      
      typedef utils::chart<active_set_type, std::allocator<active_set_type> > active_chart_type;
      typedef std::vector<active_chart_type, std::allocator<active_chart_type> > active_chart_set_type;
      
      typedef utils::chunk_vector<passive_type, 4096 / sizeof(passive_type), std::allocator<passive_type> > passive_set_type;
      typedef utils::chart<passive_set_type, std::allocator<passive_set_type> > passive_chart_type;
      
      typedef google::dense_hash_map<symbol_type, int, boost::hash<symbol_type>, std::equal_to<symbol_type> > node_map_type;
      
      typedef std::vector<unary_type, std::allocator<unary_type> > unary_set_type;
#ifdef HAVE_TR1_UNORDERED_MAP
      typedef std::tr1::unordered_map<symbol_type, unary_set_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				      std::allocator<std::pair<const symbol_type, unary_set_type> > > unary_map_type;
#else
      typedef sgi::hash_map<symbol_type, unary_set_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
			    std::allocator<std::pair<const symbol_type, unary_set_type> > > unary_map_type;
#endif
      typedef google::dense_hash_map<symbol_type, score_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > closure_set_type;
      
      ParseCKY(const symbol_type& __goal,
	       const grammar_type& __grammar,
	       const function_type& __function,
	       const bool __yield_source=false,
	       const bool __treebank=false,
	       const bool __pos_mode=false)
	: goal(__goal), grammar(__grammar), function(__function), yield_source(__yield_source), treebank(__treebank), pos_mode(__pos_mode)
      {
	node_map.set_empty_key(symbol_type());
	
	closure.set_empty_key(symbol_type());
	closure_next.set_empty_key(symbol_type());
      }
      
      template <typename Pruner>
      bool operator()(const lattice_type& lattice, label_score_chart_type& scores, const Pruner& pruner)
      {
	inside.clear();
	outside.clear();
	scores.clear();
	
	inside.resize(lattice.size() + 1);
	outside.resize(lattice.size() + 1);
	scores.resize(lattice.size() + 1);
	
	actives.clear();
	passives.clear();
	passives_unary.clear();
	
	actives.resize(grammar.size(), active_chart_type(lattice.size() + 1));
	passives.resize(lattice.size() + 1);
	passives_unary.resize(lattice.size() + 1);
	
	
	compute_inside(lattice, pruner);
	
	const bool has_goal = inside(0, lattice.size()).find(goal) != inside(0, lattice.size()).end();

	if (! has_goal) return false;

	compute_outside(lattice);
	compute_inside_outside(lattice, scores);

	return true;
      }
      
      void compute_inside_outside(const lattice_type& lattice, label_score_chart_type& scores)
      {
	//std::cerr << "inside-outside" << std::endl;

	const score_type score_sum = inside(0, lattice.size())[goal];

	// we simply enumerate chart...!
	for (size_type length = 1; length <= lattice.size(); ++ length)
	  for (size_type first = 0; first + length <= lattice.size(); ++ first) {
	    const size_type last = first + length;
	    
	    label_score_set_type&  labels_scores = scores(first, last);
	    const label_score_set_type& labels_inside = inside(first, last);
	    
	    typename label_score_set_type::const_iterator oiter_end = outside(first, last).end();
	    for (typename label_score_set_type::const_iterator oiter = outside(first, last).begin(); oiter != oiter_end; ++ oiter) {
	      typename label_score_set_type::const_iterator iiter = labels_inside.find(oiter->first);
	      if (iiter != labels_inside.end())
		labels_scores[oiter->first] = oiter->second * iiter->second / score_sum;
	    }
	  }
      }
      
      void compute_outside(const lattice_type& lattice)
      {
	//std::cerr << "outside" << std::endl;

	// traverse back passives from TOP
	
	// find goal node out of passives_unary.
	//
	// how do we traverse back this complicated structure....
	//
	
	outside(0, lattice.size())[goal] = cicada::semiring::traits<score_type>::one();
	
	for (size_type length = lattice.size(); length != 0; -- length)
	  for (size_type first = 0; first + length <= lattice.size(); ++ first) {
	    const size_type last = first + length;
	    
	    // first, enumerate unary rules
	    const passive_set_type& unaries = passives_unary(first, last);
	    
	    typename passive_set_type::const_iterator uiter_end = unaries.end();
	    for (typename passive_set_type::const_iterator uiter = unaries.begin(); uiter != uiter_end; ++ uiter) {
	      const symbol_type& head = uiter->span.label;
	      
	      typename label_score_set_type::const_iterator oiter = outside(first, last).find(head);
	      if (oiter == outside(first, last).end()) continue;

	      const score_type score_head = oiter->second;
	      
	      typename edge_set_type::const_iterator eiter_end = uiter->edges.end();
	      for (typename edge_set_type::const_iterator eiter = uiter->edges.begin(); eiter != eiter_end; ++ eiter) {
		const edge_type& edge = *eiter;
		const symbol_type& tail = edge.tails.front().label;
		
		score_type& score = outside(first, last)[tail];
		
		score = std::max(score, score_head * edge.score);
	      }
	    }
	    
	    // second, enumerate non-unary rules
	    const passive_set_type& rules = passives(first, last);
	    
	    typename passive_set_type::const_iterator riter_end = rules.end();
	    for (typename passive_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	      const symbol_type& head = riter->span.label;
	      
	      typename label_score_set_type::const_iterator oiter = outside(first, last).find(head);
	      if (oiter == outside(first, last).end()) continue;
	      
	      const score_type score_head = oiter->second;
	      
	      typename edge_set_type::const_iterator eiter_end = riter->edges.end();
	      for (typename edge_set_type::const_iterator eiter = riter->edges.begin(); eiter != eiter_end; ++ eiter) {
		const edge_type& edge = *eiter;
		
		const score_type score_edge = score_head * edge.score;
		
		typename tail_set_type::const_iterator titer_end = edge.tails.end();
		for (typename tail_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer) {
		  score_type score_outside = score_edge;
		  for (typename tail_set_type::const_iterator niter = edge.tails.begin(); niter != titer_end; ++ niter)
		    if (titer != niter)
		      score_outside *= inside(niter->first, niter->last)[niter->label];
		  
		  score_type& score = outside(titer->first, titer->last)[titer->label];
		  score = std::max(score, score_outside);
		}
	      }
	    }
	  }
      }
      
      template <typename Pruner>
      void compute_inside(const lattice_type& lattice, const Pruner& pruner)
      {
	//std::cerr << "inside: " << grammar.size() << std::endl;

	// initialize active chart...
	for (size_type table = 0; table != grammar.size(); ++ table) {
	  const transducer_type::id_type root = grammar[table].root();
	  
	  for (size_type pos = 0; pos != lattice.size(); ++ pos)
	    if (grammar[table].valid_span(pos, pos, 0))
	      actives[table](pos, pos).push_back(active_type(root));
	}
	
	// compute inside...
	for (size_type length = 1; length <= lattice.size(); ++ length)
	  for (size_type first = 0; first + length <= lattice.size(); ++ first) {
	    const size_type last = first + length;

	    // check pruning!
	    if (pruner(first, last)) continue;
	    
	    //std::cerr << "span: " << first << ".." << last << std::endl;
	    
	    node_map.clear();
	    
	    for (size_t table = 0; table != grammar.size(); ++ table) {
	      const transducer_type& transducer = grammar[table];
	      
	      // we will advance active spans, but constrained by transducer's valid span
	      if (transducer.valid_span(first, last, lattice.shortest_distance(first, last))) {
		active_set_type& cell = actives[table](first, last);
		for (size_t middle = first + 1; middle < last; ++ middle) {
		  const active_set_type&  active_arcs  = actives[table](first, middle);
		  const passive_set_type& passive_arcs = passives_unary(middle, last);
		  
		  extend_actives(transducer, active_arcs, passive_arcs, cell);
		}
		
		if (! treebank || length == 1) {
		  // then, advance by terminal(s) at lattice[last - 1];
		  
		  const active_set_type&  active_arcs  = actives[table](first, last - 1);
		  const lattice_type::arc_set_type& passive_arcs = lattice[last - 1];
		  
		  typename active_set_type::const_iterator aiter_begin = active_arcs.begin();
		  typename active_set_type::const_iterator aiter_end = active_arcs.end();
		  
		  if (aiter_begin != aiter_end) {
		    if (pos_mode) {
		      lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
		      for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) {
			const symbol_type terminal = piter->label.terminal();
			
			active_set_type& cell = actives[table](first, last - 1 + piter->distance);
			
			// handling of EPSILON rule...
			if (terminal == vocab_type::EPSILON) {
			  for (typename active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
			    cell.push_back(active_type(aiter->node, edge_type(aiter->edge.tails, aiter->edge.score + function(piter->features))));
			} else {
			  for (typename active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
			    const transducer_type::id_type node = transducer.next(aiter->node, terminal);
			    if (node == transducer.root()) continue;
			    
			    cell.push_back(active_type(node, edge_type(aiter->edge.tails, aiter->edge.score * function(piter->features))));
			  }
			}
		      }
		    } else {
		      lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
		      for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) {
			const symbol_type& terminal = piter->label;
			
			active_set_type& cell = actives[table](first, last - 1 + piter->distance);
			
			// handling of EPSILON rule...
			if (terminal == vocab_type::EPSILON) {
			  for (typename active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
			    cell.push_back(active_type(aiter->node, edge_type(aiter->edge.tails, aiter->edge.score + function(piter->features))));
			} else {
			  for (typename active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
			    const transducer_type::id_type node = transducer.next(aiter->node, terminal);
			    if (node == transducer.root()) continue;
			    
			    cell.push_back(active_type(node, edge_type(aiter->edge.tails, aiter->edge.score * function(piter->features))));
			  }
			}
		      }
		    }
		  }
		}
	      }
	      
	      // complete active items if possible... The active items may be created from child span due to the
	      // lattice structure...
	      // apply rules on actives at [first, last)
	      
	      active_set_type&  cell         = actives[table](first, last);
	      passive_set_type& passive_arcs = passives(first, last);
	      label_score_set_type& labels_inside = inside(first, last);
	      
	      typename active_set_type::const_iterator citer_end = cell.end();
	      for (typename active_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) {
		const transducer_type::rule_pair_set_type& rules = transducer.rules(citer->node);
		
		if (rules.empty()) continue;
		
		score_type score_tails = cicada::semiring::traits<score_type>::one();
		typename tail_set_type::const_iterator titer_end = citer->edge.tails.end();
		for (typename tail_set_type::const_iterator titer = citer->edge.tails.begin(); titer != titer_end; ++ titer)
		  score_tails *= inside(titer->first, titer->last)[titer->label];
		
		transducer_type::rule_pair_set_type::const_iterator riter_begin = rules.begin();
		transducer_type::rule_pair_set_type::const_iterator riter_end   = rules.end();
		
		for (transducer_type::rule_pair_set_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter) {
		  const rule_ptr_type rule = (yield_source ? riter->source : riter->target);
		  const symbol_type& lhs = rule->lhs;
		  
		  // check pruning!
		  if (pruner(first, last, lhs)) continue;

		  //std::cerr << "rule: " << *rule << std::endl;
		  
		  std::pair<typename node_map_type::iterator, bool> result = node_map.insert(std::make_pair(lhs, passive_arcs.size()));
		  if (result.second)
		    passive_arcs.push_back(span_type(first, last, lhs));
		  
		  const score_type score_edge = citer->edge.score * function(riter->features);
		  
		  passive_arcs[result.first->second].edges.push_back(edge_type(citer->edge.tails, score_edge));
		  
		  score_type& score = labels_inside[lhs];
		  score = std::max(score, score_tails * score_edge);
		}
	      }
	    }
	    
	    if (! passives(first, last).empty()) {
	      // unary rules...
	      // we will cache unary rule application by caching closure...
	      // how to avoid cycles...?
	      //
	      // compute closure from the list of possives...
	      // we assume that the rules are maximum-likely estimated, meaning negative scores..
	      // 
	      // we will simply perform max-like computation...
	      
	      // the obvious unary chain is zero-rules...
	      // simply insert a single-tail edge, with log-score of zero (or prob of 1)
	      //

	      node_map.clear();

	      passive_set_type& passive_unary = passives_unary(first, last);
	      label_score_set_type& labels_inside = inside(first, last);
	      
	      typename passive_set_type::const_iterator piter_end = passives(first, last).end();
	      for (typename passive_set_type::const_iterator piter = passives(first, last).begin(); piter != piter_end; ++ piter) {
		// child to parent...
		const unary_set_type& closure = unary_closure(piter->span.label);

		const score_type score_tail = labels_inside[piter->span.label];
		
		typename unary_set_type::const_iterator citer_end = closure.end();
		for (typename unary_set_type::const_iterator citer = closure.begin(); citer != citer_end; ++ citer) {
		  // check pruning!
		  if (pruner(first, last, citer->label)) continue;
		  
		  std::pair<typename node_map_type::iterator, bool> result = node_map.insert(std::make_pair(citer->label, passive_unary.size()));
		  if (result.second)
		    passive_unary.push_back(span_type(first, last, citer->label));
		  
		  passive_unary[result.first->second].edges.push_back(edge_type(tail_set_type(1, piter->span), citer->score));
		  
		  //std::cerr << "unary: " << citer->label << " ||| " << piter->span.label << std::endl;
		  
		  score_type& score = labels_inside[citer->label];
		  score = std::max(score, score_tail * citer->score);
		}
	      }
	    }
	    
	    // extend root with passive items at [first, last)
	    for (size_t table = 0; table != grammar.size(); ++ table) {
	      const transducer_type& transducer = grammar[table];
	      
	      if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
	      
	      const active_set_type&  active_arcs  = actives[table](first, first);
	      const passive_set_type& passive_arcs = passives_unary(first, last);
	      
	      active_set_type& cell = actives[table](first, last);
	      
	      extend_actives(transducer, active_arcs, passive_arcs, cell);
	    }
	  }
      }
      
      const unary_set_type& unary_closure(const symbol_type& child)
      {
	
	//
	// given this child state, compute closure...
	// we do not allow cycle, and keep only max-rules
	//

	typename unary_map_type::iterator uiter = unaries.find(child);
	if (uiter == unaries.end()) {
	  closure.clear();
	  closure_next.clear();
	  closure.insert(std::make_pair(child, cicada::semiring::traits<score_type>::one()));
	  
	  for (;;) {
	    bool equilibrate = true;

	    closure_next = closure;
	    
	    typename closure_set_type::const_iterator citer_end = closure.end();
	    for (typename closure_set_type::const_iterator citer = closure.begin(); citer != citer_end; ++ citer) {
	      for (size_type table = 0; table != grammar.size(); ++ table) {
		const transducer_type& transducer = grammar[table];
		
		const transducer_type::id_type node = transducer.next(transducer.root(), citer->first);
		if (node == transducer.root()) continue;
		
		const transducer_type::rule_pair_set_type& rules = transducer.rules(node);
		
		if (rules.empty()) continue;
		
		transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
		for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
		  const rule_ptr_type rule = (yield_source ? riter->source : riter->target);
		  const symbol_type& lhs = rule->lhs;
		  
		  // we assume max-like estimate grammar!
		  if (lhs == child) continue;
		  
		  const score_type score = function(riter->features) * citer->second;
		  
		  std::pair<typename closure_set_type::iterator, bool> result = closure_next.insert(std::make_pair(lhs, score));
		  if (result.second)
		    equilibrate = false;
		  else if (result.first->second < score) {
		    equilibrate = false;
		    result.first->second = score;
		  }
		}
	      }
	    }
	    
	    closure.swap(closure_next);
	    closure_next.clear();
	    
	    if (equilibrate) break;
	  }
	  
	  uiter = unaries.insert(std::make_pair(child, unary_set_type(closure.begin(), closure.end()))).first;
	}
	return uiter->second;
      }
      
      bool extend_actives(const transducer_type& transducer,
			  const active_set_type& actives, 
			  const passive_set_type& passives,
			  active_set_type& cell)
      {
	typename active_set_type::const_iterator aiter_begin = actives.begin();
	typename active_set_type::const_iterator aiter_end   = actives.end();
	
	typename passive_set_type::const_iterator piter_begin = passives.begin();
	typename passive_set_type::const_iterator piter_end   = passives.end();
	
	bool found = false;
	
	if (piter_begin != piter_end)
	  for (typename active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	    if (transducer.has_next(aiter->node)) {
	      symbol_type label;
	      transducer_type::id_type node = transducer.root();
	      
	      tail_set_type tails(aiter->edge.tails.size() + 1);
	      std::copy(aiter->edge.tails.begin(), aiter->edge.tails.end(), tails.begin());
	      
	      for (typename passive_set_type::const_iterator piter = piter_begin; piter != piter_end; ++ piter) {
		const symbol_type& non_terminal = piter->span.label;
		
		if (label != non_terminal) {
		  node = transducer.next(aiter->node, non_terminal);
		  label = non_terminal;
		}
		if (node == transducer.root()) continue;
		
		tails.back() = piter->span;
		cell.push_back(active_type(node, edge_type(tails, aiter->edge.score)));
		
		found = true;
	      }
	    }
      
	return found;
      }
      
      const symbol_type goal;
      const grammar_type& grammar;
      
      const function_type& function;
      
      const bool yield_source;
      const bool treebank;
      const bool pos_mode;

      node_map_type node_map;
      
      label_score_chart_type inside;
      label_score_chart_type outside;
      
      active_chart_set_type actives;
      passive_chart_type passives;
      passive_chart_type passives_unary;
      
      unary_map_type unaries;
      closure_set_type closure;
      closure_set_type closure_next;
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
	pos_mode(__pos_mode)
    {
      if (grammars.empty())
	throw std::runtime_error("no grammar?");
      if (thresholds.size() + 1 != grammars.size())
	throw std::runtime_error("do we have enough threshold parameters for grammars?");
    }
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      if (lattice.empty()) return;

      label_score_chart_type scores_init;
      label_score_chart_type scores;
      label_score_chart_type scores_prev;
      
      //std::cerr << "initial" << std::endl;
      {
	ParseCKY parser(goal, grammars.front(), function, yield_source, treebank, pos_mode);
	parser(lattice, scores_init, PruneNone());
      }
      
      double factor = 1.0;
      
      for (;;) {
	scores = scores_init;
	
	bool succeed = true;
	
	// corse-to-fine 
	for (size_t level = 1; level != grammars.size() - 1; ++ level) {
	  //std::cerr << "level: " << level << std::endl;

	  ParseCKY parser(goal, grammars[level], function, yield_source, treebank, pos_mode);
	  
	  scores_prev.swap(scores);
	  
	  if (level == 1)
	    succeed = parser(lattice, scores, PruneCoarse<CoarseSimple>(scores_prev, thresholds[level - 1] * factor, CoarseSimple()));
	  else
	    succeed = parser(lattice, scores, PruneCoarse<CoarseSymbol>(scores_prev, thresholds[level - 1] * factor, CoarseSymbol(level - 2)));
	  
	  if (! succeed) break;
	}
	
	if (! succeed) {
	  factor *= 0.5;
	  continue;
	}
	
	// final parsing with hypergraph construction
	ComposeCKY composer(goal, grammars.back(), yield_source, treebank, pos_mode, true);
	
	// we will fallback to simple tag!
	if (grammars.size() == 2)
	  composer(lattice, graph, PruneCoarse<CoarseSimple>(scores,
							     thresholds.back() * factor,
							     CoarseSimple()));
	else
	  composer(lattice, graph, PruneCoarse<CoarseSymbol>(scores,
							     thresholds.back() * factor,
							     CoarseSymbol(grammars.size() - 2)));
	
	break;
      }
    }
    
  private:
    const symbol_type goal;
    grammar_set_type   grammars;
    threshold_set_type thresholds;
    
    const function_type& function;

    const bool yield_source;
    const bool treebank;
    const bool pos_mode;
  };
  
  
  template <typename IteratorGrammar, typename IteratorThreshold, typename Function>
  inline
  void parse_coarse(const Symbol& goal, 
		    IteratorGrammar gfirst, IteratorGrammar glast,
		    IteratorThreshold tfirst, IteratorThreshold tlast,
		    const Function& function,
		    const Lattice& lattice,
		    HyperGraph& graph,
		    const bool yield_source=false,
		    const bool treebank=false,
		    const bool pos_mode=false)
  {
    ParseCoarse<typename Function::value_type, Function>(goal, gfirst, glast, tfirst, tlast, function, yield_source, treebank, pos_mode)(lattice, graph);
  }
  
  
};

#endif
