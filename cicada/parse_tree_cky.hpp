// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_TREE_CKY__HPP__
#define __CICADA__PARSE_TREE_CKY__HPP__ 1

#include <vector>
#include <algorithm>
#include <set>
#include <queue>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/tree_grammar.hpp>
#include <cicada/tree_transducer.hpp>
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

#include <google/dense_hash_map>
#include <google/dense_hash_set>

#include <boost/fusion/tuple.hpp>

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
  struct ParseTreeCKY
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;
    
    typedef TreeGrammar    tree_grammar_type;
    typedef TreeTransducer tree_transducer_type;
    
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    
    typedef Lattice    lattice_type;
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef attribute_set_type::attribute_type attribute_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    typedef tree_transducer_type::rule_type          tree_rule_type;
    typedef tree_transducer_type::rule_ptr_type      tree_rule_ptr_type;

    typedef Semiring semiring_type;
    typedef Semiring score_type;
    
    typedef Function function_type;
    
    ParseTreeCKY(const symbol_type& __goal,
		 const tree_grammar_type& __tree_grammar,
		 const grammar_type& __grammar,
		 const function_type& __function,
		 const int __beam_size,
		 const bool __yield_source=false,
		 const bool __unique_goal=false)
      : goal(__goal),
	tree_grammar(__tree_grammar),
	grammar(__grammar),
	function(__function),
	beam_size(__beam_size),
	yield_source(__yield_source),
	unique_goal(__unique_goal),
	attr_span_first("span-first"),
	attr_span_last("span-last"),
	attr_glue_tree(__grammar.empty() ? "" : "glue-tree")
    {
      goal_rule = rule_type::create(rule_type(vocab_type::GOAL, rule_type::symbol_set_type(1, goal.non_terminal())));
    }
    
    struct ActiveTree
    {
      tree_transducer_type::id_type             node;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type   features;
      attribute_set_type attributes;
      
      ActiveTree(const hypergraph_type::edge_type::node_set_type& __tails)
	: node(), tails(__tails), features(), attributes() {}
      ActiveTree(const tree_transducer_type::id_type& __node)
	: node(__node), tails(), features(), attributes() {}
      ActiveTree(const tree_transducer_type::id_type& __node,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node), tails(), features(__features), attributes(__attributes) {}
      ActiveTree(const tree_transducer_type::id_type& __node,
		 const hypergraph_type::edge_type::node_set_type& __tails,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node), tails(__tails), features(__features), attributes(__attributes) {}
    };
    
    struct ActiveRule
    {
      transducer_type::id_type                  node;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type   features;
      attribute_set_type attributes;
      
      ActiveRule(const hypergraph_type::edge_type::node_set_type& __tails)
	: node(), tails(__tails), features(), attributes() {}
      ActiveRule(const transducer_type::id_type& __node)
	: node(__node), tails(), features(), attributes() {}
      ActiveRule(const transducer_type::id_type& __node,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node), tails(), features(__features), attributes(__attributes) {}
      ActiveRule(const transducer_type::id_type& __node,
		 const hypergraph_type::edge_type::node_set_type& __tails,
		 const feature_set_type& __features,
		 const attribute_set_type& __attributes)
	: node(__node), tails(__tails), features(__features), attributes(__attributes) {}
    };

    typedef ActiveTree active_tree_type;
    typedef ActiveRule active_rule_type;
    
    typedef utils::chunk_vector<active_tree_type, 4096 / sizeof(active_tree_type), std::allocator<active_tree_type> > active_tree_set_type;
    typedef utils::chunk_vector<active_rule_type, 4096 / sizeof(active_rule_type), std::allocator<active_rule_type> > active_rule_set_type;
    
    typedef utils::chart<active_tree_set_type, std::allocator<active_tree_set_type> > active_tree_chart_type;
    typedef utils::chart<active_rule_set_type, std::allocator<active_rule_set_type> > active_rule_chart_type;
    
    typedef std::vector<active_tree_chart_type, std::allocator<active_tree_chart_type> > active_tree_chart_set_type;
    typedef std::vector<active_rule_chart_type, std::allocator<active_rule_chart_type> > active_rule_chart_set_type;

    struct TreeCandidate
    {
      score_type    score;
      
      symbol_type        lhs;
      tree_rule_ptr_type rule;
      
      feature_set_type   features;
      attribute_set_type attributes;
      
      TreeCandidate() : score(), lhs(), rule(), features(), attributes() {}
      TreeCandidate(const score_type& __score, const symbol_type& __lhs, const tree_rule_ptr_type& __rule, const feature_set_type& __features, const attribute_set_type& __attributes)
	: score(__score), lhs(__lhs), rule(__rule), features(__features), attributes(__attributes) {}

      void swap(TreeCandidate& x)
      {
	std::swap(score, x.score);
	lhs.swap(x.lhs);
	rule.swap(x.rule);
	features.swap(x.features);
	attributes.swap(x.attributes);
      }
      
      friend
      void swap(TreeCandidate& x, TreeCandidate& y)
      {
	x.swap(y);
      }
    };
    
    struct RuleCandidate
    {
      score_type    score;
      
      symbol_type   lhs;
      rule_ptr_type rule;
      
      feature_set_type   features;
      attribute_set_type attributes;
      
      RuleCandidate() : score(), rule(), features(), attributes() {}
      RuleCandidate(const score_type& __score, const symbol_type& __lhs, const rule_ptr_type& __rule, const feature_set_type& __features, const attribute_set_type& __attributes)
	: score(__score), lhs(__lhs), rule(__rule), features(__features), attributes(__attributes) {}

      
      void swap(RuleCandidate& x)
      {
	std::swap(score, x.score);
	lhs.swap(x.lhs);
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
    
    typedef TreeCandidate tree_candidate_type;
    typedef RuleCandidate rule_candidate_type;
    
    typedef utils::simple_vector<tree_candidate_type, std::allocator<tree_candidate_type> > tree_candidate_set_type;
    typedef utils::simple_vector<rule_candidate_type, std::allocator<rule_candidate_type> > rule_candidate_set_type;
    
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<tree_transducer_type::id_type, tree_candidate_set_type, utils::hashmurmur<size_t>, std::equal_to<tree_transducer_type::id_type>,
				    std::allocator<std::pair<const tree_transducer_type::id_type, tree_candidate_set_type> > > tree_candidate_map_type;
    typedef std::tr1::unordered_map<transducer_type::id_type, rule_candidate_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
				    std::allocator<std::pair<const transducer_type::id_type, rule_candidate_set_type> > > rule_candidate_map_type;
#else
    typedef sgi::hash_map<tree_transducer_type::id_type, tree_candidate_set_type, utils::hashmurmur<size_t>, std::equal_to<tree_transducer_type::id_type>,
			  std::allocator<std::pair<const tree_transducer_type::id_type, tree_candidate_set_type> > > tree_candidate_map_type;
    typedef sgi::hash_map<transducer_type::id_type, rule_candidate_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
			  std::allocator<std::pair<const transducer_type::id_type, rule_candidate_set_type> > > rule_candidate_map_type;
#endif
    typedef std::vector<tree_candidate_map_type, std::allocator<tree_candidate_map_type> > tree_candidate_table_type;
    typedef std::vector<rule_candidate_map_type, std::allocator<rule_candidate_map_type> > rule_candidate_table_type;

    
#ifdef HAVE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_set<const tree_candidate_type*, utils::hashmurmur<size_t>, std::equal_to<const tree_candidate_type*>,
				    std::allocator<const tree_candidate_type*> > unary_tree_set_type;
    typedef std::tr1::unordered_set<const rule_candidate_type*, utils::hashmurmur<size_t>, std::equal_to<const rule_candidate_type*>,
				    std::allocator<const rule_candidate_type*> > unary_rule_set_type;
#else
    typedef sgi::hash_set<const tree_candidate_type*, utils::hashmurmur<size_t>, std::equal_to<const tree_candidate_type*>,
			  std::allocator<const tree_candidate_type*> > unary_tree_set_type;
    typedef sgi::hash_set<const rule_candidate_type*, utils::hashmurmur<size_t>, std::equal_to<const rule_candidate_type*>,
			  std::allocator<const rule_candidate_type*> > unary_rule_set_type;
#endif
  
    typedef std::pair<symbol_type, int> symbol_level_type;
    typedef std::pair<symbol_level_type, symbol_level_type> symbol_level_pair_type;
    
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<symbol_level_pair_type, unary_tree_set_type, utils::hashmurmur<size_t>, std::equal_to<symbol_level_pair_type>,
				    std::allocator< std::pair<const symbol_level_pair_type, unary_tree_set_type> > > unary_tree_map_type;
    typedef std::tr1::unordered_map<symbol_level_pair_type, unary_rule_set_type, utils::hashmurmur<size_t>, std::equal_to<symbol_level_pair_type>,
				    std::allocator< std::pair<const symbol_level_pair_type, unary_rule_set_type> > > unary_rule_map_type;
#else
    typedef sgi::hash_map<symbol_level_pair_type, unary_tree_set_type, utils::hashmurmur<size_t>, std::equal_to<symbol_level_pair_type>,
			  std::allocator<std::pair<const symbol_level_pair_type, unary_tree_set_type> > > unary_tree_map_type;
    typedef sgi::hash_map<symbol_level_pair_type, unary_rule_set_type, utils::hashmurmur<size_t>, std::equal_to<symbol_level_pair_type>,
			  std::allocator<std::pair<const symbol_level_pair_type, unary_rule_set_type> > > unary_rule_map_type;
#endif
    
    struct Candidate
    {
      Candidate() : active_rule(0), active_tree(0), rule_first(), rule_last(), tree_first(), tree_last(), score(), level(0) {}
      Candidate(const active_rule_type* __active,
		typename rule_candidate_set_type::const_iterator __first,
		typename rule_candidate_set_type::const_iterator __last,
		const score_type& __score,
		const int __level=0)
	: active_rule(__active), active_tree(0),
	  rule_first(__first), rule_last(__last), tree_first(), tree_last(),
	  score(__score),
	  level(__level) {}
      Candidate(const active_tree_type* __active,
		typename tree_candidate_set_type::const_iterator __first,
		typename tree_candidate_set_type::const_iterator __last,
		const score_type& __score,
		const int __level=0)
	: active_rule(0), active_tree(__active),
	  rule_first(), rule_last(), tree_first(__first), tree_last(__last),
	  score(__score),
	  level(__level) {}
      
      const active_rule_type* active_rule;
      const active_tree_type* active_tree;
      
      typename rule_candidate_set_type::const_iterator rule_first;
      typename rule_candidate_set_type::const_iterator rule_last;
      typename tree_candidate_set_type::const_iterator tree_first;
      typename tree_candidate_set_type::const_iterator tree_last;
            
      score_type score;
      int level;

      bool is_rule() const { return active_rule; }
      bool is_tree() const { return active_tree; }

      score_type candidate_score() const
      {
	return score * (active_rule ? rule_first->score : tree_first->score);
      }
    };
    
    typedef Candidate candidate_type;
    typedef utils::chunk_vector<candidate_type, 1024 * 8 / sizeof(candidate_type), std::allocator<candidate_type> > candidate_set_type;
    
    struct compare_heap_type
    {
      // we use less, so that when popped from heap, we will grab "greater" in back...
      bool operator()(const candidate_type* x, const candidate_type* y) const
      {
	return x->candidate_score() < y->candidate_score();
      }
    };
    
    typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_base_type;
    typedef utils::std_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type> candidate_heap_type;
    
    typedef hypergraph_type::id_type passive_type;
    typedef std::vector<passive_type, std::allocator<passive_type> > passive_set_type;
    typedef utils::chart<passive_set_type, std::allocator<passive_set_type> > passive_chart_type;

    struct symbol_level_hash : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const symbol_level_type& x) const
      {
	return hasher_type::operator()(x.first, x.second);
      }
    };
    
    class NodeMap : public google::dense_hash_map<symbol_level_type, hypergraph_type::id_type, symbol_level_hash, std::equal_to<symbol_level_type> >
    {
    public:
      typedef google::dense_hash_map<symbol_level_type, hypergraph_type::id_type, symbol_level_hash, std::equal_to<symbol_level_type> > node_map_type;
      
    public:
      NodeMap() : node_map_type() { node_map_type::set_empty_key(symbol_level_type(symbol_type(), -1)); }
    };
    typedef NodeMap node_map_type;

    class NodeSet : public google::dense_hash_map<symbol_type, hypergraph_type::id_type, boost::hash<symbol_type>, std::equal_to<symbol_type> >
    {
    public:
      typedef google::dense_hash_map<symbol_type, hypergraph_type::id_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > node_set_type;
      
    public:
      NodeSet() : node_set_type() { node_set_type::set_empty_key(symbol_type()); }
    };
    typedef NodeSet node_set_type;
    
    typedef utils::chunk_vector<node_set_type, 4096 / sizeof(node_set_type), std::allocator<node_set_type> > node_graph_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > non_terminal_set_type;
    typedef std::vector<score_type,  std::allocator<score_type> >  score_set_type;

    typedef hypergraph_type::edge_type::node_set_type tail_set_type;
    typedef rule_type::symbol_set_type                symbol_set_type;
    typedef boost::fusion::tuple<tail_set_type, symbol_set_type, symbol_type> internal_label_type;

    struct internal_label_hash : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;

      size_t operator()(const internal_label_type& x) const
      {
	return hasher_type::operator()(boost::fusion::get<0>(x).begin(), boost::fusion::get<0>(x).end(),
				       hasher_type::operator()(boost::fusion::get<1>(x).begin(), boost::fusion::get<1>(x).end(), boost::fusion::get<2>(x).id()));
      }
    };

#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<internal_label_type, hypergraph_type::id_type, internal_label_hash, std::equal_to<internal_label_type>,
				    std::allocator<std::pair<const internal_label_type, hypergraph_type::id_type> > > internal_label_map_type;
#else
    typedef sgi::hash_map<internal_label_type, hypergraph_type::id_type, internal_label_hash, std::equal_to<internal_label_type>,
			  std::allocator<std::pair<const internal_label_type, hypergraph_type::id_type> > > internal_label_map_type;
#endif

    
    struct less_non_terminal
    {
      less_non_terminal(const non_terminal_set_type& __non_terminals) : non_terminals(__non_terminals) {}
      
      bool operator()(const hypergraph_type::id_type& x, const hypergraph_type::id_type& y) const
      {
	return non_terminals[x] < non_terminals[y] || (non_terminals[x] == non_terminals[y] && x < y);
      }
      
      const non_terminal_set_type& non_terminals;
    };

    struct VerifyNone
    {
      template <typename Transducer>
      bool operator()(const Transducer& transducer, const size_t first, const size_t last, const size_t distance) const
      {
	return true;
      }
    };
    
    struct VerifySpan
    {
      template <typename Transducer>
      bool operator()(const Transducer& transducer, const size_t first, const size_t last, const size_t distance) const
      {
	return transducer.valid_span(first, last, distance);
      }
    };
    
    void operator()(const lattice_type& lattice,
		    hypergraph_type& graph)
    {
      graph.clear();
      
      if (lattice.empty())
	return;
      
      // initialize internal structure...
      actives_rule.clear();
      actives_tree.clear();
      passives.clear();
      
      node_graph_rule.clear();
      node_graph_tree.clear();
      non_terminals.clear();
      scores.clear();

      label_map.clear();

      actives_tree.reserve(tree_grammar.size());
      actives_rule.reserve(grammar.size());
      
      actives_rule.resize(grammar.size(), active_rule_chart_type(lattice.size() + 1));
      actives_tree.resize(tree_grammar.size(), active_tree_chart_type(lattice.size() + 1));
      
      passives.resize(lattice.size() + 1);
    
      rule_tables.clear();
      rule_tables.reserve(grammar.size());
      rule_tables.resize(grammar.size());
      
      tree_tables.clear();
      tree_tables.reserve(tree_grammar.size());
      tree_tables.resize(tree_grammar.size());
    
      // initialize active chart

      for (size_t table = 0; table != tree_grammar.size(); ++ table) {
	const tree_transducer_type::id_type root = tree_grammar[table].root();
	
	for (size_t pos = 0; pos != lattice.size(); ++ pos)
	  actives_tree[table](pos, pos).push_back(active_tree_type(root));
      }
      
      for (size_t table = 0; table != grammar.size(); ++ table) {
	const transducer_type::id_type root = grammar[table].root();
	
	for (size_t pos = 0; pos != lattice.size(); ++ pos)
	  if (grammar[table].valid_span(pos, pos, 0))
	    actives_rule[table](pos, pos).push_back(active_rule_type(root));
      }
      
      for (size_t length = 1; length <= lattice.size(); ++ length)
	for (size_t first = 0; first + length <= lattice.size(); ++ first) {
	  const size_t last = first + length;
	  	  
	  //std::cerr << "span: " << first << ".." << last << " distance: " << lattice.shortest_distance(first, last) << std::endl;
	  
	  extend_actives(first, last, lattice, tree_grammar, actives_tree, passives, VerifyNone());
	  extend_actives(first, last, lattice, grammar,      actives_rule, passives, VerifySpan());
	  
	  // complete active items if possible... The active items may be created from child span due to the
	  // lattice structure...
	  // apply rules on actives at [first, last)
	  
	  //
	  // we will try apply rules, queue them as "candidate"
	  //
	  // when candidate is popped, create graph
	  // create new candidate with unary rule, but if it is already created, ignore!
	  // 
	  
	  actives_rule_unary.clear();
	  actives_tree_unary.clear();
	  
	  unary_rule_map.clear();
	  unary_tree_map.clear();
	  
	  node_map.clear();
	  candidates.clear();
	  heap.clear();
	  
	  for (size_t table = 0; table != tree_grammar.size(); ++ table) {
	    active_tree_set_type&  cell = actives_tree[table](first, last);

	    typename active_tree_set_type::const_iterator citer_end = cell.end();
	    for (typename active_tree_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) {
	      const tree_candidate_set_type& rules = candidate_trees(table, citer->node);
	      
	      if (rules.empty()) continue;
	      
	      score_type score_antecedent = function(citer->features);
	      
	      hypergraph_type::edge_type::node_set_type::const_iterator titer_end = citer->tails.end();
	      for (hypergraph_type::edge_type::node_set_type::const_iterator titer = citer->tails.begin(); titer != titer_end; ++ titer)
		score_antecedent *= scores[*titer];
	      
	      candidates.push_back(candidate_type(&(*citer), rules.begin(), rules.end(), score_antecedent));
	      heap.push(&candidates.back());
	    }
	  }
	  
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    active_rule_set_type&  cell = actives_rule[table](first, last);
	    
	    typename active_rule_set_type::const_iterator citer_end = cell.end();
	    for (typename active_rule_set_type::const_iterator citer = cell.begin(); citer != citer_end; ++ citer) {
	      const rule_candidate_set_type& rules = candidate_rules(table, citer->node);
	      
	      if (rules.empty()) continue;
	      
	      score_type score_antecedent = function(citer->features);
	      
	      hypergraph_type::edge_type::node_set_type::const_iterator titer_end = citer->tails.end();
	      for (hypergraph_type::edge_type::node_set_type::const_iterator titer = citer->tails.begin(); titer != titer_end; ++ titer)
		score_antecedent *= scores[*titer];
	      
	      candidates.push_back(candidate_type(&(*citer), rules.begin(), rules.end(), score_antecedent));
	      heap.push(&candidates.back());
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
	    
	    const score_type score = item->candidate_score();
	    
	    std::pair<hypergraph_type::id_type, bool> node_passive;
	    
	    if (item->is_tree()) {
	      const active_tree_type& active = *(item->active_tree);
	      const tree_candidate_type& rule = *(item->tree_first);
	      
	      if (item->level > 0) {
		const symbol_type label_prev = non_terminals[active.tails.front()];
		const symbol_type label_next = rule.lhs;
		
		unary_tree_set_type& unaries = unary_tree_map[std::make_pair(std::make_pair(label_prev, item->level - 1), std::make_pair(label_next, item->level))];
		
		if (! unaries.insert(&(*(item->tree_first))).second) {
		  typename node_map_type::const_iterator niter = node_map.find(std::make_pair(label_next, item->level));
		  if (niter == node_map.end())
		    throw std::runtime_error("no node-map?");
		  
		  node_passive.first = niter->second;
		  node_passive.second = score > scores[niter->second];
		  
		  scores[niter->second] = std::max(scores[niter->second], score);
		} else
		  node_passive = apply_rule(score,
					    rule.lhs,
					    rule.rule,
					    active.features + rule.features,
					    active.attributes + rule.attributes,
					    active.tails,
					    passive_arcs,
					    graph,
					    first,
					    last,
					    item->level);
	      } else
		node_passive = apply_rule(score,
					  rule.lhs,
					  rule.rule,
					  active.features + rule.features,
					  active.attributes + rule.attributes,
					  active.tails,
					  passive_arcs,
					  graph,
					  first,
					  last,
					  item->level);

	      // next queue!
	      ++ const_cast<candidate_type*>(item)->tree_first;
	      if (item->tree_first != item->tree_last)
		heap.push(item);
	      
	    } else {
	      const active_rule_type& active = *(item->active_rule);
	      const rule_candidate_type& rule = *(item->rule_first);
	      
	      if (item->level > 0) {
		const symbol_type label_prev = non_terminals[active.tails.front()];
		const symbol_type label_next = rule.lhs;
		
		unary_rule_set_type& unaries = unary_rule_map[std::make_pair(std::make_pair(label_prev, item->level - 1), std::make_pair(label_next, item->level))];
		if (! unaries.insert(&(*(item->rule_first))).second) {
		  typename node_map_type::const_iterator niter = node_map.find(std::make_pair(label_next, item->level));
		  if (niter == node_map.end())
		    throw std::runtime_error("no node-map?");
		  
		  node_passive.first = niter->second;
		  node_passive.second = score > scores[niter->second];
		  
		  scores[niter->second] = std::max(scores[niter->second], score);
		} else
		  node_passive = apply_rule(score,
					    rule.lhs,
					    rule.rule,
					    active.features + rule.features,
					    active.attributes + rule.attributes,
					    active.tails,
					    passive_arcs,
					    graph,
					    first,
					    last,
					    item->level);
	      } else
		node_passive = apply_rule(score,
					  rule.lhs,
					  rule.rule,
					  active.features + rule.features,
					  active.attributes + rule.attributes,
					  active.tails,
					  passive_arcs,
					  graph,
					  first,
					  last,
					  item->level);
	      
	      // next queue!
	      ++ const_cast<candidate_type*>(item)->rule_first;
	      if (item->rule_first != item->rule_last)
		heap.push(item);
	    }
	    
	    if (! node_passive.second) continue;
	    
	    const symbol_type& non_terminal = non_terminals[node_passive.first];
	    const score_type score_antecedent = scores[node_passive.first];
	    
	    // extend by unary rule...
	    for (size_t table = 0; table != tree_grammar.size(); ++ table) {
	      const tree_transducer_type& transducer = tree_grammar[table];
	      
	      const tree_transducer_type::id_type node = transducer.next(transducer.root(), non_terminal);
	      if (node == transducer.root()) continue;
	      
	      const tree_candidate_set_type& rules = candidate_trees(table, node);
	      
	      if (rules.empty()) continue;
	      
	      actives_tree_unary.push_back(active_tree_type(hypergraph_type::edge_type::node_set_type(1, node_passive.first)));
	      candidates.push_back(candidate_type(&actives_tree_unary.back(), rules.begin(), rules.end(), score_antecedent, item->level + 1));
	      heap.push(&candidates.back());
	    }

	    for (size_t table = 0; table != grammar.size(); ++ table) {
	      const transducer_type& transducer = grammar[table];
	      
	      if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
	      
	      const transducer_type::id_type node = transducer.next(transducer.root(), non_terminal);
	      if (node == transducer.root()) continue;
	      
	      const rule_candidate_set_type& rules = candidate_rules(table, node);
	      
	      if (rules.empty()) continue;
	      
	      actives_rule_unary.push_back(active_rule_type(hypergraph_type::edge_type::node_set_type(1, node_passive.first)));
	      candidates.push_back(candidate_type(&actives_rule_unary.back(), rules.begin(), rules.end(), score_antecedent, item->level + 1));
	      heap.push(&candidates.back());
	    }
	  }
	  
	  // sort passives at passives(first, last) wrt non-terminal label in non_terminals
	  {
	    passive_set_type& passive_arcs = passives(first, last);
	    
	    passive_set_type(passive_arcs).swap(passive_arcs);
	    std::sort(passive_arcs.begin(), passive_arcs.end(), less_non_terminal(non_terminals));
	  }

	  //std::cerr << "span: " << first << ".." << last << " passives: " << passives(first, last).size() << std::endl;
	  
	  // extend root with passive items at [first, last)
	  for (size_t table = 0; table != tree_grammar.size(); ++ table) {
	    const tree_transducer_type& transducer = tree_grammar[table];
	    
	    const active_tree_set_type& active_arcs  = actives_tree[table](first, first);
	    const passive_set_type&     passive_arcs = passives(first, last);
	    
	    active_tree_set_type& cell = actives_tree[table](first, last);
	    
	    extend_actives(transducer, active_arcs, passive_arcs, cell);
	  }
	  
	  for (size_t table = 0; table != grammar.size(); ++ table) {
	    const transducer_type& transducer = grammar[table];
	    
	    if (! transducer.valid_span(first, last, lattice.shortest_distance(first, last))) continue;
	    
	    const active_rule_set_type&  active_arcs  = actives_rule[table](first, first);
	    const passive_set_type&      passive_arcs = passives(first, last);
	    
	    active_rule_set_type& cell = actives_rule[table](first, last);
	    
	    extend_actives(transducer, active_arcs, passive_arcs, cell);
	  }
	}
      
      //
      // final patch work..
      //
      // we will connect node_graph_rule into node_graph_tree
      //
      for (size_t node_id = 0; node_id != node_graph_tree.size(); ++ node_id) {
	const node_set_type& node_set_tree = node_graph_tree[node_id];
	const node_set_type& node_set_rule = node_graph_rule[node_id];
	
	if (node_set_tree.empty() || node_set_rule.empty()) continue;

	typename node_set_type::const_iterator citer_begin = node_set_rule.begin();
	typename node_set_type::const_iterator citer_end   = node_set_rule.end();

	typename node_set_type::const_iterator piter_end = node_set_tree.end();
	for (typename node_set_type::const_iterator piter = node_set_tree.begin(); piter != piter_end; ++ piter)
	  for (typename node_set_type::const_iterator citer = citer_begin; citer != citer_end; ++ citer) {
	    hypergraph_type::edge_type& edge = graph.add_edge(&(citer->second), &(citer->second) + 1);
	    
	    edge.rule = rule_type::create(rule_type(piter->first, &(citer->first), &(citer->first) + 1));
	    edge.attributes[attr_glue_tree] = attribute_set_type::int_type(1);
	    
	    graph.connect_edge(edge.id, piter->second);
	  }
      }
      
      // final goal assignment...
      if (unique_goal) {
	passive_set_type& passive_arcs = passives(0, lattice.size());
	for (size_t p = 0; p != passive_arcs.size(); ++ p) {
	  const node_set_type& node_set = node_graph_tree[passive_arcs[p]];
	  
	  typename node_set_type::const_iterator giter = node_set.find(goal);
	  if (giter != node_set.end()) {
	    if (graph.is_valid())
	      throw std::runtime_error("multiple goal? " + boost::lexical_cast<std::string>(graph.goal) + " " + boost::lexical_cast<std::string>(passive_arcs[p]));
	    
	    graph.goal = giter->second;
	  }
	}
	
      } else {
	passive_set_type& passive_arcs = passives(0, lattice.size());
	for (size_t p = 0; p != passive_arcs.size(); ++ p) {
	  const node_set_type& node_set = node_graph_tree[passive_arcs[p]];
	  
	  typename node_set_type::const_iterator giter = node_set.find(goal);
	  if (giter == node_set.end()) continue;
	  
	  hypergraph_type::edge_type& edge = graph.add_edge(&(giter->second), &(giter->second) + 1);
	  edge.rule = goal_rule;
	  edge.attributes[attr_span_first] = attribute_set_type::int_type(0);
	  edge.attributes[attr_span_last]  = attribute_set_type::int_type(lattice.size());
	  
	  if (! graph.is_valid())
	    graph.goal = graph.add_node().id;
	  
	  graph.connect_edge(edge.id, graph.goal);
	}
      }
      
      // we will sort to remove unreachable nodes......
      if (graph.is_valid())
	graph.topologically_sort();
    }

  private:
    
    std::pair<hypergraph_type::id_type, bool> apply_rule(const score_type& score,
							 const symbol_type& lhs,
							 const rule_ptr_type& rule,
							 const feature_set_type& features,
							 const attribute_set_type& attributes,
							 const hypergraph_type::edge_type::node_set_type& frontier,
							 passive_set_type& passives,
							 hypergraph_type& graph,
							 const int lattice_first,
							 const int lattice_last,
							 const int level = 0)
    {
       //
      // we need to transform the source-frontier into target-frontier!
      //
      
      hypergraph_type::edge_type::node_set_type tails(frontier.size());
      if (! frontier.empty()) {
	int non_terminal_pos = 0;
	rule_type::symbol_set_type::const_iterator riter_end = rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator riter = rule->rhs.begin(); riter != riter_end; ++ riter)
	  if (riter->is_non_terminal()) {
	    const int __non_terminal_index = riter->non_terminal_index();
	    const int non_terminal_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	    ++ non_terminal_pos;
	    
	    node_set_type& node_set = node_graph_rule[frontier[non_terminal_index]];
	    std::pair<typename node_set_type::iterator, bool> result = node_set.insert(std::make_pair(riter->non_terminal(), 0));
	    if (result.second)
	      result.first->second = graph.add_node().id;
	    
	    tails[non_terminal_index] = result.first->second;
	  }
      }
      
      hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
      edge.rule = rule;
      edge.features = features;
      edge.attributes = attributes;
      
      // assign metadata...
      edge.attributes[attr_span_first] = attribute_set_type::int_type(lattice_first);
      edge.attributes[attr_span_last]  = attribute_set_type::int_type(lattice_last);
      
      bool unary_next = false;

      std::pair<typename node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(lhs, level), 0));
      if (result.second) {
	result.first->second = node_graph_tree.size();
	
	non_terminals.push_back(lhs);
	passives.push_back(node_graph_tree.size());
	scores.push_back(score);
	
	node_graph_tree.resize(node_graph_tree.size() + 1);
	node_graph_rule.resize(node_graph_rule.size() + 1);
	
	unary_next = true;
      } else
	unary_next = score > scores[result.first->second];
      
      scores[result.first->second] = std::max(scores[result.first->second], score);
      
      // projected lhs
      std::pair<typename node_set_type::iterator, bool> result_mapped = node_graph_rule[result.first->second].insert(std::make_pair(rule->lhs, 0));
      if (result_mapped.second)
	result_mapped.first->second = graph.add_node().id;
      
      graph.connect_edge(edge.id, result_mapped.first->second);
      
      return std::make_pair(result.first->second, unary_next);
    }


    std::pair<hypergraph_type::id_type, bool> apply_rule(const score_type& score,
							 const symbol_type& lhs,
							 const tree_rule_ptr_type& rule,
							 const feature_set_type& features,
							 const attribute_set_type& attributes,
							 const hypergraph_type::edge_type::node_set_type& frontier,
							 passive_set_type& passives,
							 hypergraph_type& graph,
							 const int lattice_first,
							 const int lattice_last,
							 const int level = 0)
    {
      bool unary_next = false;
      
      // source lhs
      std::pair<typename node_map_type::iterator, bool> result = node_map.insert(std::make_pair(std::make_pair(lhs, level), 0));
      if (result.second) {
	result.first->second = node_graph_tree.size();
	
	non_terminals.push_back(lhs);
	passives.push_back(node_graph_tree.size());
	scores.push_back(score);
	
	node_graph_tree.resize(node_graph_tree.size() + 1);
	node_graph_rule.resize(node_graph_rule.size() + 1);

	unary_next = true;
      } else
	unary_next = score > scores[result.first->second];
      
      scores[result.first->second] = std::max(scores[result.first->second], score);
      
      // projected lhs
      std::pair<typename node_set_type::iterator, bool> result_mapped = node_graph_tree[result.first->second].insert(std::make_pair(rule->label, 0));
      if (result_mapped.second)
	result_mapped.first->second = graph.add_node().id;
      
      const hypergraph_type::id_type root_id = result_mapped.first->second;
      
      int non_terminal_pos = 0;
      
      const hypergraph_type::id_type edge_id = construct_graph(*rule, root_id, frontier, graph, non_terminal_pos);
      
      graph.edges[edge_id].features   = features;
      graph.edges[edge_id].attributes = attributes;
      
      // assign metadata only for the root edge...?????
      graph.edges[edge_id].attributes[attr_span_first] = attribute_set_type::int_type(lattice_first);
      graph.edges[edge_id].attributes[attr_span_last]  = attribute_set_type::int_type(lattice_last);
      
      return std::make_pair(result.first->second, unary_next);
    }

    hypergraph_type::id_type construct_graph(const tree_rule_type& rule,
					     hypergraph_type::id_type root,
					     const hypergraph_type::edge_type::node_set_type& frontiers,
					     hypergraph_type& graph,
					     int& non_terminal_pos)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > rhs_type;
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails_type;
      
      rhs_type rhs;
      tails_type tails;
      
      tree_rule_type::const_iterator aiter_end = rule.end();
      for (tree_rule_type::const_iterator aiter = rule.begin(); aiter != aiter_end; ++ aiter)
	if (aiter->label.is_non_terminal()) {
	  if (aiter->antecedents.empty()) {
	    const int __non_terminal_index = aiter->label.non_terminal_index();
	    const int non_terminal_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	    ++ non_terminal_pos;
	    
	    if (non_terminal_index >= static_cast<int>(frontiers.size()))
	      throw std::runtime_error("non-terminal index exceeds frontier size");
	    
	    node_set_type& node_set = node_graph_tree[frontiers[non_terminal_index]];
	    std::pair<typename node_set_type::iterator, bool> result = node_set.insert(std::make_pair(aiter->label.non_terminal(), 0));
	    if (result.second)
	      result.first->second = graph.add_node().id;
	    
	    // transform into frontier of the translational forest
	    tails.push_back(result.first->second);
	  } else
	    tails.push_back(0);
	  
	  rhs.push_back(aiter->label.non_terminal());
	} else
	  rhs.push_back(aiter->label);

      // if we have antecedents traverse and construct
      tails_type::iterator titer = tails.begin();
      for (tree_rule_type::const_iterator aiter = rule.begin(); aiter != aiter_end; ++ aiter)
	if (aiter->label.is_non_terminal()) {
	  if (! aiter->antecedents.empty()) {
	    const hypergraph_type::id_type edge_id = construct_graph(*aiter, hypergraph_type::invalid, frontiers, graph, non_terminal_pos);
	    *titer = graph.edges[edge_id].head;
	  }
	  ++ titer;
	}

      hypergraph_type::id_type edge_id;
      
      if (root == hypergraph_type::invalid) {
#if 1
	// we will share internal nodes
	
	std::pair<typename internal_label_map_type::iterator, bool> result = label_map.insert(std::make_pair(internal_label_type(tail_set_type(tails.begin(), tails.end()),
																 symbol_set_type(rhs.begin(), rhs.end()),
																 rule.label), 0));
	if (result.second) {
	  root = graph.add_node().id;
	  
	  edge_id = graph.add_edge(tails.begin(), tails.end()).id;
	  graph.edges[edge_id].rule = rule_type::create(rule_type(rule.label, rhs.begin(), rhs.end()));
	  graph.connect_edge(edge_id, root);
	  
	  result.first->second = edge_id;
	} else
	  edge_id = result.first->second;
#endif
#if 0
	root = graph.add_node().id;
	edge_id = graph.add_edge(tails.begin(), tails.end()).id;
	graph.edges[edge_id].rule = rule_type::create(rule_type(rule.label, rhs.begin(), rhs.end()));
	graph.connect_edge(edge_id, root);
#endif
      } else {
	edge_id = graph.add_edge(tails.begin(), tails.end()).id;
	graph.edges[edge_id].rule = rule_type::create(rule_type(rule.label, rhs.begin(), rhs.end()));
	graph.connect_edge(edge_id, root);
      }
      
      return edge_id;
    }
    
    template <typename Transducers, typename Actives, typename Verify>
    void extend_actives(const size_t first,
			const size_t last,
			const lattice_type& lattice,
			const Transducers& transducers,
			Actives& actives,
			const passive_chart_type& passives,
			Verify verify)
    {
      typedef typename Actives::value_type::value_type active_set_type;
      typedef typename Transducers::transducer_type transducer_type;

      for (size_t table = 0; table != transducers.size(); ++ table) {
	const transducer_type& transducer = transducers[table];
	
	if (! verify(transducer, first, last, lattice.shortest_distance(first, last))) continue;
	
	active_set_type& cell = actives[table](first, last);
	
	for (size_t middle = first + 1; middle < last; ++ middle)
	  extend_actives(transducer, actives[table](first, middle), passives(middle, last), cell);
	
	// then, advance by terminal(s) at lattice[last - 1];
	const active_set_type&            active_arcs  = actives[table](first, last - 1);
	const lattice_type::arc_set_type& passive_arcs = lattice[last - 1];
	
	typename active_set_type::const_iterator aiter_begin = active_arcs.begin();
	typename active_set_type::const_iterator aiter_end   = active_arcs.end();
	
	if (aiter_begin == aiter_end) continue;
	
	lattice_type::arc_set_type::const_iterator piter_end = passive_arcs.end();
	for (lattice_type::arc_set_type::const_iterator piter = passive_arcs.begin(); piter != piter_end; ++ piter) {
	  const symbol_type& terminal = piter->label;
	  
	  active_set_type& cell = actives[table](first, last - 1 + piter->distance);
	  
	  // handling of EPSILON tree...
	  if (terminal == vocab_type::EPSILON) {
	    for (typename active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	      cell.push_back(typename active_set_type::value_type(aiter->node, aiter->tails, aiter->features + piter->features, aiter->attributes));
	  } else {
	    for (typename active_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
	      const typename transducer_type::id_type node = transducer.next(aiter->node, terminal);
	      if (node == transducer.root()) continue;
	      
	      cell.push_back(typename active_set_type::value_type(node, aiter->tails, aiter->features + piter->features, aiter->attributes));
	    }
	  }
	}
      }
    }
    
    template <typename Transducer, typename Actives>
    bool extend_actives(const Transducer& transducer,
			const Actives& actives, 
			const passive_set_type& passives,
			Actives& cell)
    {
      typename Actives::const_iterator aiter_begin = actives.begin();
      typename Actives::const_iterator aiter_end   = actives.end();
      
      passive_set_type::const_iterator piter_begin = passives.begin();
      passive_set_type::const_iterator piter_end   = passives.end();
      
      bool found = false;
      
      if (piter_begin != piter_end)
	for (typename Actives::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  if (transducer.has_next(aiter->node)) {
	    symbol_type label;
	    typename Transducer::id_type node = transducer.root();
	    
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
	      cell.push_back(typename Actives::value_type(node, tails, aiter->features, aiter->attributes));
	      
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
    
    const rule_candidate_set_type& candidate_rules(const size_type& table, const transducer_type::id_type& node)
    {
      typename rule_candidate_map_type::iterator riter = rule_tables[table].find(node);
      if (riter == rule_tables[table].end()) {
	const transducer_type::rule_pair_set_type& rules = grammar[table].rules(node);
	
	riter = rule_tables[table].insert(std::make_pair(node, rule_candidate_set_type(rules.size()))).first;
	
	typename rule_candidate_set_type::iterator citer = riter->second.begin();
	transducer_type::rule_pair_set_type::const_iterator iter_end   = rules.end();
	for (transducer_type::rule_pair_set_type::const_iterator iter = rules.begin(); iter != iter_end; ++ iter, ++ citer)
	  *citer = rule_candidate_type(function(iter->features),
				       iter->source->lhs,
				       yield_source ? iter->source : iter->target,
				       iter->features,
				       iter->attributes);
	  
	std::sort(riter->second.begin(), riter->second.end(), greater_score<rule_candidate_type>());
      }
      
      return riter->second;
    }

    const tree_candidate_set_type& candidate_trees(const size_type& table, const tree_transducer_type::id_type& node)
    {
      typename tree_candidate_map_type::iterator riter = tree_tables[table].find(node);
      if (riter == tree_tables[table].end()) {
	const tree_transducer_type::rule_pair_set_type& rules = tree_grammar[table].rules(node);
	
	riter = tree_tables[table].insert(std::make_pair(node, tree_candidate_set_type(rules.size()))).first;
	
	
	typename tree_candidate_set_type::iterator citer = riter->second.begin();
	tree_transducer_type::rule_pair_set_type::const_iterator iter_end   = rules.end();
	for (tree_transducer_type::rule_pair_set_type::const_iterator iter = rules.begin(); iter != iter_end; ++ iter, ++ citer)
	  *citer = tree_candidate_type(function(iter->features),
				       iter->source->label,
				       yield_source ? iter->source : iter->target,
				       iter->features,
				       iter->attributes);
	
	std::sort(riter->second.begin(), riter->second.end(), greater_score<tree_candidate_type>());
      }
      
      return riter->second;
    }
    
  private:
    const symbol_type goal;
    const tree_grammar_type& tree_grammar;
    const grammar_type& grammar;
    
    const function_type& function;
    const int beam_size;

    const bool yield_source;
    const bool unique_goal;
    const attribute_type attr_span_first;
    const attribute_type attr_span_last;
    const attribute_type attr_glue_tree;
    
    rule_ptr_type goal_rule;

    active_rule_chart_set_type actives_rule;
    active_tree_chart_set_type actives_tree;
    passive_chart_type         passives;
    active_rule_set_type       actives_rule_unary;
    active_tree_set_type       actives_tree_unary;

    candidate_set_type    candidates;
    candidate_heap_type   heap;

    rule_candidate_table_type rule_tables;
    tree_candidate_table_type tree_tables;
    
    unary_rule_map_type   unary_rule_map;
    unary_tree_map_type   unary_tree_map;
    node_map_type         node_map;
    node_graph_type       node_graph_tree;
    node_graph_type       node_graph_rule;
    non_terminal_set_type non_terminals;
    score_set_type        scores;

    internal_label_map_type label_map;
  };
  
  template <typename Function>
  inline
  void parse_tree_cky(const Symbol& goal, const TreeGrammar& tree_grammar, const Grammar& grammar, const Function& function, const Lattice& lattice, HyperGraph& graph, const int size, const bool yield_source=false, const bool unique_goal=false)
  {
    ParseTreeCKY<typename Function::value_type, Function>(goal, tree_grammar, grammar, function, size, yield_source, unique_goal)(lattice, graph);
  }
  
};

#endif
