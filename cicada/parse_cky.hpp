// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_CKY__HPP__
#define __CICADA__PARSE_CKY__HPP__ 1

#include <vector>
#include <algorithm>
#include <set>
#include <queue>

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
#include <utils/simple_vector.hpp>
#include <utils/small_vector.hpp>
#include <utils/mulvector2.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

namespace cicada
{
  // semiring and function to compute semiring from a feature vector
  
  // CKY + beam search
  // Basic idea is to collect all the active items wrt passive items,
  // queue them during rule-application, and grab the most promissing rules...
  // do cube-pruning with passive items?
  // TODO: how to handle unary chains?
  // hash-based hypergraph structure, then, finally convert into actual graph???
  // 

  template <typename Semiring, typename Function>
  struct ParseCKY
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
    
    ParseCKY(const symbol_type& __goal,
	     const grammar_type& __grammar,
	     const function_type& __function,
	     const int __beam_size,
	     const bool __yield_source=false,
	     const bool __treebank=false,
	     const bool __pos_mode=false,
	     const bool __unique_goal=false)
      : goal(__goal),
	grammar(__grammar),
	function(__function),
	beam_size(__beam_size),
	yield_source(__yield_source),
	treebank(__treebank),
	pos_mode(__pos_mode),
	unique_goal(__unique_goal),
	attr_span_first("span-first"),
	attr_span_last("span-last")
    {
      goal_rule = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, goal.non_terminal())));
      
      node_map.set_empty_key(symbol_level_type());
    }
    
    struct Active
    {
      Active(const transducer_type::id_type& __node,
	     const hypergraph_type::edge_type::node_set_type& __tails,
	     const feature_set_type& __features,
	     const attribute_set_type& __attributes)
	: node(__node),
	  tails(__tails),
	  features(__features),
	  attributes(__attributes) {}
      Active(const transducer_type::id_type& __node,
	     const feature_set_type& __features,
	     const attribute_set_type& __attributes)
	: node(__node),
	  tails(),
	  features(__features),
	  attributes(__attributes) {}
      Active(const transducer_type::id_type& __node)
	: node(__node),
	  tails(),
	  features(),
	  attributes() {}
      
      Active()
	: node(),
	  tails(),
	  features(),
	  attributes() {}
      
      transducer_type::id_type                  node;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type                          features;
      attribute_set_type                        attributes;
    };
    
    typedef Active active_type;
    typedef utils::chunk_vector<active_type, 4096 / sizeof(active_type), std::allocator<active_type> > active_set_type;

    typedef utils::chart<active_set_type, std::allocator<active_set_type> > active_chart_type;
    typedef std::vector<active_chart_type, std::allocator<active_chart_type> > active_chart_set_type;

    struct RuleCandidate
    {
      score_type    score;
      
      rule_ptr_type rule;
      
      feature_set_type   features;
      attribute_set_type attributes;
      
      RuleCandidate() : score(), rule(), features(), attributes() {}
      RuleCandidate(const score_type& __score, const rule_ptr_type& __rule, const feature_set_type& __features, const attribute_set_type& __attributes)
	: score(__score), rule(__rule), features(__features), attributes(__attributes) {}

      void swap(RuleCandidate& x)
      {
	std::swap(score, x.score);
	rule.swap(x.rule);
	features.swap(x.features);
	attributes.swap(x.attributes);
      }
      
      friend
      void swap(RuleCandidate& x, RuleCandidate& y)
      {
	x.swap(y);
      }
    };
    typedef RuleCandidate rule_candidate_type;
    typedef utils::simple_vector<rule_candidate_type, std::allocator<rule_candidate_type> > rule_candidate_set_type;
    
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<transducer_type::id_type, rule_candidate_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
				    std::allocator<std::pair<const transducer_type::id_type, rule_candidate_set_type> > > rule_candidate_map_type;
#else
    typedef sgi::hash_map<transducer_type::id_type, rule_candidate_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
			  std::allocator<std::pair<const transducer_type::id_type, rule_candidate_set_type> > > rule_candidate_map_type;
#endif
    typedef std::vector<rule_candidate_map_type, std::allocator<rule_candidate_map_type> > rule_candidate_table_type;
    
    
#ifdef HAVE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_set<const rule_candidate_type*, utils::hashmurmur<size_t>, std::equal_to<const rule_candidate_type*>,
				    std::allocator<const rule_candidate_type*> > unary_rule_set_type;
#else
    typedef sgi::hash_set<const rule_candidate_type*, utils::hashmurmur<size_t>, std::equal_to<const rule_candidate_type*>,
			  std::allocator<const rule_candidate_type*> > unary_rule_set_type;
    
#endif
    
    typedef std::pair<symbol_type, int> symbol_level_type;
    
    struct symbol_level_hash : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const symbol_level_type& x) const
      {
	return hasher_type::operator()(x.first, x.second);
      }
    };
    
    typedef google::dense_hash_map<symbol_level_type, hypergraph_type::id_type, symbol_level_hash, std::equal_to<symbol_level_type> > node_map_type;
    
    typedef std::pair<symbol_level_type, symbol_level_type> symbol_level_pair_type;

#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<symbol_level_pair_type, unary_rule_set_type, utils::hashmurmur<size_t>, std::equal_to<symbol_level_pair_type>,
				    std::allocator< std::pair<const symbol_level_pair_type, unary_rule_set_type> > > unary_rule_map_type;
#else
    typedef sgi::hash_map<symbol_level_pair_type, unary_rule_set_type, utils::hashmurmur<size_t>, std::equal_to<symbol_level_pair_type>,
			  std::allocator<std::pair<const symbol_level_pair_type, unary_rule_set_type> > > unary_rule_map_type;

#endif
    
    typedef utils::small_vector<int, std::allocator<int> > index_set_type;
    
    struct Candidate
    {
      Candidate() : active(), j(), first(), iter(), last(), score(), level(0) {}
      
      const active_type* active;
      
      index_set_type j;
      
      typename rule_candidate_set_type::const_iterator first;
      typename rule_candidate_set_type::const_iterator iter;
      typename rule_candidate_set_type::const_iterator last;

      score_type score;
      int level;
    };
    typedef Candidate candidate_type;
    typedef utils::chunk_vector<candidate_type, 1024 * 8 / sizeof(candidate_type), std::allocator<candidate_type> > candidate_set_type;
    
    struct compare_heap_type
    {
      // we use less, so that when popped from heap, we will grab "greater" in back...
      bool operator()(const candidate_type* x, const candidate_type* y) const
      {
	//return x->score * x->iter->score < y->score * y->iter->score;
	return x->score < y->score;
      }
    };
    
    typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_base_type;
    typedef utils::std_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type> candidate_heap_type;
    
    typedef hypergraph_type::id_type passive_type;
    typedef std::vector<passive_type, std::allocator<passive_type> > passive_set_type;
    typedef utils::chart<passive_set_type, std::allocator<passive_set_type> > passive_chart_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
    typedef std::vector<score_type,  std::allocator<score_type> >  score_set_type;
    
    struct less_non_terminal
    {
      less_non_terminal(const non_terminal_set_type& __non_terminals,
			const score_set_type& __scores)
	: non_terminals(__non_terminals),
	  scores(__scores) {}
      
      bool operator()(const hypergraph_type::id_type& x, const hypergraph_type::id_type& y) const
      {
	return (non_terminals[x] < non_terminals[y]
		|| (non_terminals[x] == non_terminals[y]
		    && (scores[x] > scores[y]
			|| (scores[x] == scores[y]
			    && x < y))));
      }
      
      const non_terminal_set_type& non_terminals;
      const score_set_type&        scores;
    };
    
    typedef std::vector<passive_type, std::allocator<passive_type> > derivation_set_type;
    typedef utils::mulvector2<passive_type, std::allocator<passive_type> > derivation_map_type;
    
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
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      operator()(lattice, graph, PruneNone());
    }
    
    template <typename Pruner>
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph,
		    const Pruner& pruner)
    {
      graph.clear();
      
      if (lattice.empty())
	return;
      
      // initialize internal structure...
      actives.clear();
      passives.clear();
      non_terminals.clear();
      scores.clear();
      derivation_map.clear();

      actives.reserve(grammar.size());
      passives.reserve(lattice.size() + 1);
    
      actives.resize(grammar.size(), active_chart_type(lattice.size() + 1));
      passives.resize(lattice.size() + 1);

      rule_tables.clear();
      rule_tables.reserve(grammar.size());
      rule_tables.resize(grammar.size());

      derivation_set_type derivations;
    
      // initialize active chart
      for (size_t table = 0; table != grammar.size(); ++ table) {
	const transducer_type::id_type root = grammar[table].root();
	
	for (size_t pos = 0; pos != lattice.size(); ++ pos)
	  if (grammar[table].valid_span(pos, pos, 0))
	    actives[table](pos, pos).push_back(active_type(root));
      }
      
      for (size_t length = 1; length <= lattice.size(); ++ length)
	for (size_t first = 0; first + length <= lattice.size(); ++ first) {
	  const size_t last = first + length;

	  if (pruner(first, last)) continue;
	  	  
	  //std::cerr << "span: " << first << ".." << last << " distance: " << lattice.shortest_distance(first, last) << std::endl;
	  
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	    
	    // we will advance active spans, but constrained by transducer's valid span
	    if (transducer.valid_span(first, last, lattice.shortest_distance(first, last))) {
	      // advance dots....
	      
	      // first, extend active items...
	      active_set_type& cell = actives[table](first, last);
	      for (size_t middle = first + 1; middle < last; ++ middle) {
		const active_set_type&  active_arcs  = actives[table](first, middle);
		const passive_set_type& passive_arcs = passives(middle, last);
		
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
			  cell.push_back(active_type(aiter->node, aiter->tails, aiter->features + piter->features, aiter->attributes));
		      } else {
			for (typename active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
			  const transducer_type::id_type node = transducer.next(aiter->node, terminal);
			  if (node == transducer.root()) continue;
			  
			  cell.push_back(active_type(node, aiter->tails, aiter->features + piter->features, aiter->attributes));
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
			  cell.push_back(active_type(aiter->node, aiter->tails, aiter->features + piter->features, aiter->attributes));
		      } else {
			for (typename active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
			  const transducer_type::id_type node = transducer.next(aiter->node, terminal);
			  if (node == transducer.root()) continue;
			  
			  cell.push_back(active_type(node, aiter->tails, aiter->features + piter->features, aiter->attributes));
			}
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
	  
	  //
	  // we will try apply rules, queue them as "candidate"
	  //
	  // when candidate is popped, create graph
	  // create new candidate with unary rule, but if it is already created, ignore!
	  // 
	  
	  actives_unary.clear();
	  unary_map.clear();
	  node_map.clear();
	  candidates.clear();
	  heap.clear();
	  
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    active_set_type&  cell = actives[table](first, last);
	    
	    typename active_set_type::const_iterator citer_end = cell.end();
	    for (typename active_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) {
	      const rule_candidate_set_type& rules = cands(table, citer->node);
	      
	      if (rules.empty()) continue;
	      
	      score_type score_antecedent = function(citer->features);
	      
	      hypergraph_type::edge_type::node_set_type::const_iterator titer_end = citer->tails.end();
	      for (hypergraph_type::edge_type::node_set_type::const_iterator titer = citer->tails.begin(); titer != titer_end; ++ titer)
		score_antecedent *= scores[derivation_map[*titer].front()];
	      
	      candidates.push_back(candidate_type());
	      candidate_type& cand = candidates.back();
	      
	      cand.active = &(*citer);
	      cand.j = index_set_type(citer->tails.size(), 0);
	      
	      cand.first = rules.begin();
	      cand.iter  = rules.begin();
	      cand.last  = rules.end();
	      
	      cand.score = score_antecedent * cand.iter->score;
	      cand.level = 0;
	      
	      heap.push(&cand);
	    }
	  }

	  derivations.clear();
	  
	  for (int num_pop = 0; ! heap.empty() && num_pop != beam_size; /**/) {
	    // pop-best...
	    const candidate_type* item = heap.top();
	    heap.pop();
	    
	    // add into graph...
	    
	    //
	    // we will always expand into unary rules, in order to find out better unary chain!
	    //
	    
	    // check unary rule, and see if this edge is already inserted!
	    
	    const active_type& active = *(item->active);
	    const rule_candidate_type& rule = *(item->iter);
	    const score_type score = item->score;
	    
	    if (pruner(first, last, rule.rule->lhs)) {
	      // next queue!
	      push_succ(item);
	      continue;
	    }
	    
	    // we will increment here!
	    ++ num_pop;
	    
	    std::pair<hypergraph_type::id_type, bool> node_passive(hypergraph_type::invalid, false);
	    
	    if (item->level > 0) {
	      // check an edge consisting of:
	      //
	      // non_terminals[item->edge.tails.front()] (item->level - 1)
	      // item->rule->lhs (item->level)
	      // with item->table, item->node, item->pos
	      //
	      
	      // if already inserted, check node-map and update scores!
	      
	      // for level > 0, we do not use j!
	      const symbol_type label_prev = non_terminals[active.tails.front()];
	      const symbol_type label_next = rule.rule->lhs;
	      
	      unary_rule_set_type& unaries = unary_map[std::make_pair(std::make_pair(label_prev, item->level - 1), std::make_pair(label_next, item->level))];
	      
	      if (! unaries.insert(&(*(item->iter))).second) {
		typename node_map_type::const_iterator niter = node_map.find(std::make_pair(label_next, item->level));
		if (niter == node_map.end())
		  throw std::runtime_error("no node-map?");
		
		node_passive.first = niter->second;
		node_passive.second = score > scores[niter->second];
		
		scores[niter->second] = std::max(scores[niter->second], score);
	      } else
		node_passive = apply_rule(score, rule.rule, active.features + rule.features, active.attributes + rule.attributes,
					  active.tails.begin(), active.tails.end(), derivations, graph,
					  first, last, utils::bithack::branch(unique_goal && rule.rule->lhs == goal, 0, item->level));
	    } else {
	      // transform acive.tails into tails! if j is not empty!

	      if (item->j.empty())
		node_passive = apply_rule(score, rule.rule, active.features + rule.features, active.attributes + rule.attributes,
					  active.tails.begin(), active.tails.end(), derivations, graph,
					  first, last, item->level);
	      else {
		// transform active.tails into actual tails
		hypergraph_type::edge_type::node_set_type tails(active.tails);
		for (size_t i = 0; i != tails.size(); ++ i)
		  tails[i] = derivation_map[active.tails[i]][item->j[i]];
		
		node_passive = apply_rule(score, rule.rule, active.features + rule.features, active.attributes + rule.attributes,
					  tails.begin(), tails.end(), derivations, graph,
					  first, last, item->level);
	      }
	    }
	    
	    // next queue!
	    push_succ(item);
	    
	    // apply unary rule
	    //
	    // we will apply unary rule, when: new node is inserted or better score was found...
	    //

	    if (! node_passive.second) continue;
	    
	    const symbol_type& non_terminal = non_terminals[node_passive.first];
	    const score_type score_antecedent = scores[node_passive.first];
	    
	    for (size_t table = 0; table != grammar.size(); ++ table) {
	      const transducer_type& transducer = grammar[table];
	      
	      if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
	      
	      const transducer_type::id_type node = transducer.next(transducer.root(), non_terminal);
	      if (node == transducer.root()) continue;
	      
	      const rule_candidate_set_type& rules = cands(table, node);
	      
	      if (rules.empty()) continue;
	      
	      //std::cerr << "unary rule: " << non_terminal << " size: " << rules.size() << std::endl;
	      
	      actives_unary.push_back(active_type());
	      actives_unary.back().tails = hypergraph_type::edge_type::node_set_type(1, node_passive.first);
	      
	      candidates.push_back(candidate_type());
	      candidate_type& cand = candidates.back();
	      
	      cand.active = &(actives_unary.back());
	      // cand.j is empty!
	      cand.first = rules.begin();
	      cand.iter  = rules.begin();
	      cand.last  = rules.end();
	      
	      cand.score = score_antecedent * cand.iter->score;
	      cand.level = item->level + 1;
	      
	      heap.push(&cand);
	    }
	  }
	  
	  if (! derivations.empty()) {
	    if (derivations.size() == 1)
	      passives(first, last).push_back(derivation_map.push_back(derivations.begin(), derivations.end()));
	    else {
	      std::sort(derivations.begin(), derivations.end(), less_non_terminal(non_terminals, scores));
	      
	      //std::cerr << "derivations: " << derivations.size() << std::endl;
	      
	      // construct passives!
	      passive_set_type& passive_arcs = passives(first, last);
	      
	      size_t i_first = 0;
	      for (size_t i = 1; i != derivations.size(); ++ i)
		if (non_terminals[derivations[i_first]] != non_terminals[derivations[i]]) {
		  passive_arcs.push_back(derivation_map.push_back(derivations.begin() + i_first, derivations.begin() + i));

		  //std::cerr << "\trange: [" << i_first << ", " << i << ")" << std::endl;
		  
		  i_first = i;
		}
	      
	      if (i_first != derivations.size()) {
		passive_arcs.push_back(derivation_map.push_back(derivations.begin() + i_first, derivations.end()));
		
		//std::cerr << "\trange: [" << i_first << ", " << derivations.size() << ")" << std::endl;
	      }
	    }
	    
	    // extend root with passive items at [first, last)
	    for (size_t table = 0; table != grammar.size(); ++ table) {
	      const transducer_type& transducer = grammar[table];
	      
	      if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
	      
	      const active_set_type&  active_arcs  = actives[table](first, first);
	      const passive_set_type& passive_arcs = passives(first, last);
	      
	      active_set_type& cell = actives[table](first, last);
	      
	      extend_actives(transducer, active_arcs, passive_arcs, cell);
	    }
	  }
	  
	  //std::cerr << "span: " << first << ".." << last << " passives: " << passives(first, last).size() << std::endl;
	}
      
      // finally, collect all the parsed rules, and proceed to [goal] rule...
      // passive arcs will not be updated!
      
      // revise this!
      
      if (unique_goal) {
	passive_set_type& passive_arcs = passives(0, lattice.size());
	for (size_t p = 0; p != passive_arcs.size(); ++ p)
	  if (non_terminals[derivation_map[passive_arcs[p]].front()] == goal) {
	    if (graph.is_valid())
	      throw std::runtime_error("multiple goal? " + boost::lexical_cast<std::string>(graph.goal) + " " + boost::lexical_cast<std::string>(derivation_map[passive_arcs[p]].front()));
	    
	    if (derivation_map[passive_arcs[p]].size() > 1)
	      throw std::runtime_error("multiple goal???");
	    
	    graph.goal = derivation_map[passive_arcs[p]].front();
	  }
      } else {
	passive_set_type& passive_arcs = passives(0, lattice.size());
	for (size_t p = 0; p != passive_arcs.size(); ++ p)
	  if (non_terminals[derivation_map[passive_arcs[p]].front()] == goal) {
	    //std::cerr << "goal node: " << passive_arcs[p] << std::endl;
	    
	    typename derivation_map_type::const_reference ref = derivation_map[passive_arcs[p]];
	    
	    for (size_t i = 0; i != ref.size(); ++ i) {
	      const passive_type& id = ref[i];
	      
	      hypergraph_type::edge_type& edge = graph.add_edge(&id, (&id) + 1);
	      edge.rule = goal_rule;
	      
	      edge.attributes[attr_span_first] = attribute_set_type::int_type(0);
	      edge.attributes[attr_span_last]  = attribute_set_type::int_type(lattice.size());
	      
	      if (! graph.is_valid()) {
		graph.goal = graph.add_node().id;
		non_terminals.push_back(goal_rule->lhs);
	      }
	      
	      graph.connect_edge(edge.id, graph.goal);
	    }
	  }
      }
      
      // we will sort to remove unreachable nodes......
      if (graph.is_valid())
	graph.topologically_sort();
    }

  private:

    void push_succ(const candidate_type* item)
    {
      if (item->j.empty() || item->iter != item->first) {
	if (item->iter + 1 != item->last) {
	  const_cast<candidate_type*>(item)->score /= item->iter->score;
	  ++ const_cast<candidate_type*>(item)->iter;
	  const_cast<candidate_type*>(item)->score *= item->iter->score;
	  heap.push(item);
	}
      } else {
	if (item->iter + 1 != item->last) {
	  // make new candidate... we do not have to make an adjustment for scores...
	  candidates.push_back(*item);
	  candidate_type& cand = candidates.back();
	  cand.score /= cand.iter->score;
	  ++ cand.iter;
	  cand.score *= cand.iter->score;
	  heap.push(&cand);
	}
	
	for (size_type i = 0; i != item->j.size(); ++ i) {
	  typename derivation_map_type::const_reference ref = derivation_map[item->active->tails[i]];
	  
	  if (item->j[i] + 1 < static_cast<int>(ref.size())) {
	    // make new candidate
	    candidates.push_back(*item);
	    candidate_type& cand = candidates.back();
	    
	    ++ cand.j[i];
	    
	    // we need to adjust scores!
	    cand.score *= scores[ref[cand.j[i]]] / scores[ref[cand.j[i] - 1]];
	    
	    heap.push(&cand);
	  }
	  
	  if (item->j[i]) break;
	}
	
	// item is simply discarded!
      }
    }
    
    template <typename Iterator>
    std::pair<hypergraph_type::id_type, bool> apply_rule(const score_type& score,
							 const rule_ptr_type& rule,
							 const feature_set_type& features,
							 const attribute_set_type& attributes,
							 Iterator first,
							 Iterator last,
							 derivation_set_type& derivations,
							 hypergraph_type& graph,
							 const int lattice_first,
							 const int lattice_last,
							 const int level = 0)
    {
      //std::cerr << "rule: " << *rule << std::endl;

      hypergraph_type::edge_type& edge = graph.add_edge(first, last);
      edge.rule = rule;
      edge.features = features;
      edge.attributes = attributes;
      
      // assign metadata...
      edge.attributes[attr_span_first] = attribute_set_type::int_type(lattice_first);
      edge.attributes[attr_span_last]  = attribute_set_type::int_type(lattice_last);

      bool unary_next = false;

      std::pair<typename node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(rule->lhs, level), 0));
      if (result.second) {
	hypergraph_type::node_type& node = graph.add_node();
	
	derivations.push_back(node.id);
	non_terminals.push_back(rule->lhs);
	scores.push_back(score);
	  
	result.first->second = node.id;
	  
	unary_next = true;
      } else
	unary_next = score > scores[result.first->second];
      
      scores[result.first->second] = std::max(scores[result.first->second], score);
      
      graph.connect_edge(edge.id, result.first->second);
      
      return std::make_pair(result.first->second, unary_next);
    }
    
    bool extend_actives(const transducer_type& transducer,
			const active_set_type& actives, 
			const passive_set_type& passives,
			active_set_type& cell)
    {
      typename active_set_type::const_iterator aiter_begin = actives.begin();
      typename active_set_type::const_iterator aiter_end   = actives.end();
      
      passive_set_type::const_iterator piter_begin = passives.begin();
      passive_set_type::const_iterator piter_end   = passives.end();

      bool found = false;
      
      if (piter_begin != piter_end)
	for (typename active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  if (transducer.has_next(aiter->node)) {
	    symbol_type label;
	    transducer_type::id_type node = transducer.root();
	    
	    hypergraph_type::edge_type::node_set_type tails(aiter->tails.size() + 1);
	    std::copy(aiter->tails.begin(), aiter->tails.end(), tails.begin());
	    
	    for (passive_set_type::const_iterator piter = piter_begin; piter != piter_end; ++ piter) {
	      const symbol_type& non_terminal = non_terminals[derivation_map[*piter].front()];
	      
	      if (label != non_terminal) {
		node = transducer.next(aiter->node, non_terminal);
		label = non_terminal;
	      }
	      if (node == transducer.root()) continue;
	      
	      tails.back() = *piter;
	      cell.push_back(active_type(node, tails, aiter->features, aiter->attributes));
	      
	      found = true;
	    }
	  }
      
      return found;
    }
    
    template <typename Tp>
    struct greater_score
    {
      bool operator()(const Tp& x, const Tp& y) const
      {
	return x.score > y.score;
      }
    };

    const rule_candidate_set_type& cands(const size_type& table, const transducer_type::id_type& node)
    {
      typename rule_candidate_map_type::iterator riter = rule_tables[table].find(node);
      if (riter == rule_tables[table].end()) {
	const transducer_type::rule_pair_set_type& rules = grammar[table].rules(node);
	
	riter = rule_tables[table].insert(std::make_pair(node, rule_candidate_set_type(rules.size()))).first;
	
	typename rule_candidate_set_type::iterator citer = riter->second.begin();
	transducer_type::rule_pair_set_type::const_iterator iter_end   = rules.end();
	for (transducer_type::rule_pair_set_type::const_iterator iter = rules.begin(); iter != iter_end; ++ iter, ++ citer)
	  *citer = rule_candidate_type(function(iter->features),
				       yield_source ? iter->source : iter->target,
				       iter->features,
				       iter->attributes);
	
	std::sort(riter->second.begin(), riter->second.end(), greater_score<rule_candidate_type>());
      }
      
      return riter->second;
    }
    
  private:
    const symbol_type goal;
    const grammar_type& grammar;
    
    const function_type& function;
    const int beam_size;

    const bool yield_source;
    const bool treebank;
    const bool pos_mode;
    const bool unique_goal;
    const attribute_type attr_span_first;
    const attribute_type attr_span_last;
    
    rule_ptr_type goal_rule;

    active_chart_set_type  actives;
    passive_chart_type     passives;
    active_set_type        actives_unary;
    derivation_map_type    derivation_map;

    unary_rule_map_type   unary_map;
    node_map_type         node_map;
    candidate_set_type    candidates;
    candidate_heap_type   heap;

    rule_candidate_table_type rule_tables;
    
    non_terminal_set_type non_terminals;
    score_set_type        scores;
  };
  
  template <typename Function>
  inline
  void parse_cky(const Symbol& goal, const Grammar& grammar, const Function& function, const Lattice& lattice, HyperGraph& graph, const int size, const bool yield_source=false, const bool treebank=false, const bool pos_mode=false, const bool unique_goal=false)
  {
    ParseCKY<typename Function::value_type, Function>(goal, grammar, function, size, yield_source, treebank, pos_mode, unique_goal)(lattice, graph);
  }
  
};

#endif
