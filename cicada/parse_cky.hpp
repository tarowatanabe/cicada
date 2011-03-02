// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_CKY__HPP__
#define __CICADA__PARSE_CKY__HPP__ 1

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
	     const bool __unique_goal=false)
      : goal(__goal),
	grammar(__grammar),
	function(__function),
	beam_size(__beam_size),
	yield_source(__yield_source),
	treebank(__treebank),
	unique_goal(__unique_goal),
	attr_span_first("span-first"),
	attr_span_last("span-last")
    {
      goal_rule = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, goal.non_terminal())));
      
      node_map.set_empty_key(symbol_level_type());
    }
    
    struct ActiveItem
    {
      ActiveItem(const transducer_type::id_type& __node,
		 const hypergraph_type::edge_type::node_set_type& __tails,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node),
	  tails(__tails),
	  features(__features),
	  attributes(__attributes) {}
      ActiveItem(const transducer_type::id_type& __node,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node),
	  tails(),
	  features(__features),
	  attributes(__attributes) {}
      ActiveItem(const transducer_type::id_type& __node)
	: node(__node),
	  tails(),
	  features(),
	  attributes() {}
      
      ActiveItem()
	: node(),
	  tails(),
	  features(),
	  attributes() {}
      
      transducer_type::id_type                  node;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type                          features;
      attribute_set_type                        attributes;
    };
    
    typedef ActiveItem active_type;
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
    };
    typedef RuleCandidate rule_candidate_type;
    typedef utils::chunk_vector<rule_candidate_type, 4096 / sizeof(rule_candidate_type), std::allocator<rule_candidate_type> > rule_candidate_set_type;

    typedef std::vector<const rule_candidate_type*, std::allocator<const rule_candidate_type*> > rule_candidate_ptr_set_type;
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<transducer_type::id_type, rule_candidate_ptr_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
				    std::allocator<std::pair<const transducer_type::id_type, rule_candidate_ptr_set_type> > > rule_candidate_map_type;
#else
    typedef sgi::hash_map<transducer_type::id_type, rule_candidate_ptr_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
			  std::allocator<std::pair<const transducer_type::id_type, rule_candidate_ptr_set_type> > > rule_candidate_map_type;
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
				    std::allocator<std::pair<const symbol_level_pair_type, unary_rule_set_type> > > unary_rule_map_type;
#else
    typedef sgi::hash_map<symbol_level_pair_type, unary_rule_set_type, utils::hashmurmur<size_t>, std::equal_to<symbol_level_pair_type>,
			  std::allocator<std::pair<const symbol_level_pair_type, unary_rule_set_type> > > unary_rule_map_type;
  
#endif
    
    struct Candidate
    {
      Candidate() : active(), first(), last(), score(), level(0) {}
      
      const active_type* active;
      
      typename rule_candidate_ptr_set_type::const_iterator first;
      typename rule_candidate_ptr_set_type::const_iterator last;
      
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
	return x->score * (*(x->first))->score < y->score * (*(y->first))->score;
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
      less_non_terminal(const non_terminal_set_type& __non_terminals) : non_terminals(__non_terminals) {}

      bool operator()(const hypergraph_type::id_type& x, const hypergraph_type::id_type& y) const
      {
	return non_terminals[x] < non_terminals[y] || (non_terminals[x] == non_terminals[y] && x < y);
      }
      
      const non_terminal_set_type& non_terminals;
    };
    
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      if (lattice.empty())
	return;
      
      // initialize internal structure...
      actives.clear();
      passives.clear();
      non_terminals.clear();
      scores.clear();
    
      actives.resize(grammar.size(), active_chart_type(lattice.size() + 1));
      passives.resize(lattice.size() + 1);
    
      rule_candidates.clear();
      rule_tables.clear();
      rule_tables.reserve(grammar.size());
      rule_tables.resize(grammar.size());
    
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
	    const transducer_type& transducer = grammar[table];
	    
	    active_set_type&  cell         = actives[table](first, last);
	    
	    typename active_set_type::const_iterator citer_end = cell.end();
	    for (typename active_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) {
	      
	      const rule_candidate_ptr_set_type& rules = cands(table, citer->node);
	      
	      if (rules.empty()) continue;
	      
	      score_type score_antecedent = semiring::traits<score_type>::one();
	      
	      hypergraph_type::edge_type::node_set_type::const_iterator titer_end = citer->tails.end();
	      for (hypergraph_type::edge_type::node_set_type::const_iterator titer = citer->tails.begin(); titer != titer_end; ++ titer)
		score_antecedent *= scores[*titer];
	      
	      candidates.push_back(candidate_type());
	      candidate_type& cand = candidates.back();
	      
	      cand.active = &(*citer);
	      cand.first = rules.begin();
	      cand.last = rules.end();
	      
	      cand.score = score_antecedent;
	      cand.level = 0;
	      
	      heap.push(&cand);
	    }
	  }

	  passive_set_type& passive_arcs = passives(first, last);
	  
	  for (int num_pop = 0; ! heap.empty() && num_pop != beam_size; ++ num_pop) {
	    // pop-best...
	    const candidate_type* item = heap.top();
	    heap.pop();
	    
	    // add into graph...
	    
	    //
	    // we will always expand into unary rules, in order to find out better unary chain!
	    //
	    
	    // check unary rule, and see if this edge is already inserted!
	    
	    const active_type& active = *(item->active);
	    const rule_candidate_type& rule = *(*(item->first));
	    const score_type score = item->score * rule.score;
	    
	    const int level_next = utils::bithack::branch(unique_goal && rule.rule->lhs == goal, 0, item->level + 1);
	    
	    std::pair<hypergraph_type::id_type, bool> node_passive;
	    
	    if (item->level > 0) {
	      // check an edge consisting of:
	      //
	      // non_terminals[item->edge.tails.front()] (item->level - 1)
	      // item->rule->lhs (item->level)
	      // with item->table, item->node, item->pos
	      //
	      
	      // if already inserted, check node-map and update scores!
	      
	      const symbol_type label_prev = non_terminals[active.tails.front()];
	      const symbol_type label_next = rule.rule->lhs;
	      
	      unary_rule_set_type& unaries = unary_map[std::make_pair(std::make_pair(label_prev, item->level - 1), std::make_pair(label_next, item->level))];
	      
	      if (unaries.find(&rule) != unaries.end()) {
		typename node_map_type::const_iterator niter = node_map.find(std::make_pair(label_next, item->level));
		if (niter == node_map.end())
		  throw std::runtime_error("no node-map?");
		
		node_passive.first = niter->second;
		node_passive.second = score > scores[niter->second];
		
		scores[niter->second] = std::max(scores[niter->second], score);
	      } else {
		node_passive = apply_rule(score, rule.rule, active.features + rule.features, active.attributes + rule.attributes,
					  active.tails.begin(), active.tails.end(), passive_arcs, graph,
					  first, last, item->level);
		
		unaries.insert(&rule);
	      }
	    } else
	      node_passive = apply_rule(score, rule.rule, active.features + rule.features, active.attributes + rule.attributes,
					active.tails.begin(), active.tails.end(), passive_arcs, graph,
					first, last, item->level);
	    
	    
	    // next queue!
	    ++ const_cast<candidate_type*>(item)->first;
	    if (item->first != item->last)
	      heap.push(item);
	    
	    // apply unary rule
	    //
	    // we will apply unary rule, when: new node is inserted or better score was found...
	    //

	    if (! node_passive.second) continue;
	    
	    const symbol_type& non_terminal = non_terminals[node_passive.first];
	    
	    for (size_t table = 0; table != grammar.size(); ++ table) {
	      const transducer_type& transducer = grammar[table];
	      
	      if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
	      
	      const transducer_type::id_type node = transducer.next(transducer.root(), non_terminal);
	      if (node == transducer.root()) continue;
	      
	      const rule_candidate_ptr_set_type& rules = cands(table, node);
	      
	      if (rules.empty()) continue;
	      
	      const score_type score_antecedent = scores[node_passive.first];
	      
	      actives_unary.push_back(active_type());
	      actives_unary.back().tails = hypergraph_type::edge_type::node_set_type(1, node_passive.first);
	      
	      candidates.push_back(candidate_type());
	      candidate_type& cand = candidates.back();
	      
	      cand.active = &(actives_unary.back());
	      cand.first = rules.begin();
	      cand.last = rules.end();
	      
	      cand.score = score_antecedent;
	      cand.level = level_next;
	      
	      heap.push(&cand);
	    }
	  }
	  
	  // sort passives at passives(first, last) wrt non-terminal label in non_terminals
	  std::sort(passives(first, last).begin(), passives(first, last).end(), less_non_terminal(non_terminals));

	  //std::cerr << "span: " << first << ".." << last << " passives: " << passives(first, last).size() << std::endl;
	  
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
      
      // finally, collect all the parsed rules, and proceed to [goal] rule...
      // passive arcs will not be updated!
      
      // we will clear node map so that we will always create new node..
      node_map.clear();

      if (unique_goal) {
	passive_set_type& passive_arcs = passives(0, lattice.size());
	for (size_t p = 0; p != passive_arcs.size(); ++ p)
	  if (non_terminals[passive_arcs[p]] == goal) {
	    if (graph.is_valid())
	      throw std::runtime_error("multiple goal? " + boost::lexical_cast<std::string>(graph.goal) + " " + boost::lexical_cast<std::string>(passive_arcs[p]));
	    
	    graph.goal = passive_arcs[p];
	  }
      } else {
	passive_set_type& passive_arcs = passives(0, lattice.size());
	for (size_t p = 0; p != passive_arcs.size(); ++ p)
	  if (non_terminals[passive_arcs[p]] == goal) {
	    //std::cerr << "goal node: " << passive_arcs[p] << std::endl;
	    
	    apply_rule(score_type(), goal_rule, feature_set_type(), attribute_set_type(), &(passive_arcs[p]), (&passive_arcs[p]) + 1, passive_arcs, graph,
		       0, lattice.size(),
		       0, true);
	  }
      }
      
      // we will sort to remove unreachable nodes......
      graph.topologically_sort();
    }

  private:
    
    template <typename Iterator>
    std::pair<hypergraph_type::id_type, bool> apply_rule(const score_type& score,
							 const rule_ptr_type& rule,
							 const feature_set_type& features,
							 const attribute_set_type& attributes,
							 Iterator first,
							 Iterator last,
							 passive_set_type& passives,
							 hypergraph_type& graph,
							 const int lattice_first,
							 const int lattice_last,
							 const int level = 0,
							 const bool is_goal = false)
    {
      //std::cerr << "rule: " << *rule << std::endl;

      hypergraph_type::edge_type& edge = graph.add_edge(first, last);
      edge.rule = rule;
      edge.features = features;
      edge.attributes = attributes;
      
      // assign metadata...
      edge.attributes[attr_span_first] = attribute_set_type::int_type(lattice_first);
      edge.attributes[attr_span_last]  = attribute_set_type::int_type(lattice_last);

      if (is_goal) {
	if (! graph.is_valid()) {
	  graph.goal = graph.add_node().id;
	  non_terminals.push_back(rule->lhs);
	}
	
	graph.connect_edge(edge.id, graph.goal);
	
	return std::make_pair(graph.goal, false);
      } else {
	bool unary_next = false;

	std::pair<typename node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(rule->lhs, level), 0));
	if (result.second) {
	  hypergraph_type::node_type& node = graph.add_node();
	  
	  non_terminals.push_back(rule->lhs);
	  passives.push_back(node.id);
	  scores.push_back(score);
	  
	  result.first->second = node.id;
	  
	  unary_next = true;
	} else
	  unary_next = score > scores[result.first->second];
	
	scores[result.first->second] = std::max(scores[result.first->second], score);
	
	graph.connect_edge(edge.id, result.first->second);
	
	return std::make_pair(result.first->second, unary_next);
      }
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
	      const symbol_type& non_terminal = non_terminals[*piter];
	      
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
    struct greater_ptr_score
    {
      bool operator()(const Tp* x, const Tp* y) const
      {
	return x->score > y->score;
      }
    };

    const rule_candidate_ptr_set_type& cands(const size_type& table, const transducer_type::id_type& node)
    {
      typename rule_candidate_map_type::iterator riter = rule_tables[table].find(node);
      if (riter == rule_tables[table].end()) {
	riter = rule_tables[table].insert(std::make_pair(node, rule_candidate_ptr_set_type())).first;
	
	const transducer_type::rule_pair_set_type& rules = grammar[table].rules(node);
	
	riter->second.reserve(rules.size());
	
	transducer_type::rule_pair_set_type::const_iterator iter_begin = rules.begin();
	transducer_type::rule_pair_set_type::const_iterator iter_end   = rules.end();
	for (transducer_type::rule_pair_set_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter) {
	  rule_candidates.push_back(rule_candidate_type(function(iter->features),
							yield_source ? iter->source : iter->target,
							iter->features,
							iter->attributes));
	  
	  riter->second.push_back(&(rule_candidates.back()));
	}
	
	std::sort(riter->second.begin(), riter->second.end(), greater_ptr_score<rule_candidate_type>());
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
    const bool unique_goal;
    const attribute_type attr_span_first;
    const attribute_type attr_span_last;
    
    rule_ptr_type goal_rule;

    active_chart_set_type  actives;
    passive_chart_type     passives;
    active_set_type        actives_unary;

    unary_rule_map_type   unary_map;
    node_map_type         node_map;
    candidate_set_type    candidates;
    candidate_heap_type   heap;

    rule_candidate_set_type   rule_candidates;
    rule_candidate_table_type rule_tables;
    
    non_terminal_set_type non_terminals;
    score_set_type        scores;
  };
  
  template <typename Function>
  inline
  void parse_cky(const Symbol& goal, const Grammar& grammar, const Function& function, const Lattice& lattice, HyperGraph& graph, const int size, const bool yield_source=false, const bool treebank=false, const bool unique_goal=false)
  {
    ParseCKY<typename Function::value_type, Function>(goal, grammar, function, size, yield_source, treebank, unique_goal)(lattice, graph);
  }
  
};

#endif
