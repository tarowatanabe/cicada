// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_TREE__HPP__
#define __CICADA__PARSE_TREE__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>
#include <queue>

#include <cicada/symbol.hpp>
#include <cicada/sentence.hpp>
#include <cicada/vocab.hpp>
#include <cicada/tree_grammar.hpp>
#include <cicada/tree_transducer.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/b_heap.hpp>
#include <utils/std_heap.hpp>
#include <utils/bithack.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

#include <boost/fusion/tuple.hpp>


//
// hypergraph-to-hypergraph transduction by
//
//@InProceedings{zhang-EtAl:2009:EMNLP1,
//    author    = {Zhang, Hui  and  Zhang, Min  and  Li, Haizhou  and  Tan, Chew Lim},
//    title     = {Fast Translation Rule Matching for Syntax-based Statistical Machine Translation},
//    booktitle = {Proceedings of the 2009 Conference on Empirical Methods in Natural Language Processing},
//    month     = {August},
//    year      = {2009},
//    address   = {Singapore},
//    publisher = {Association for Computational Linguistics},
//    pages     = {1037--1045},
//    url       = {http://www.aclweb.org/anthology/D/D09/D09-1108}
//}
//
// The terminologies used in the above paper is somewhat confusing:
//   input-forest, input-hypergraph, encoded hyper-path etc.
//
// The algorithm computes:
//  for each node, try matching with tree fragments
//  when matched, the match is represented by a set of hypergraph-node-id of input-hypergraph
//  and a set of matched rules.
//  We can uncover output hypergraph by enumerating matched rules.
//  Book-keep the matched input-hypergraph in a chart so that we can construct translational packed forest.
// 

//
// beam parse variant of tree-composition
//

namespace cicada
{
  
  template <typename Semiring, typename Function>
  struct ParseTree
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;
    
    typedef Sentence sentence_type;
    typedef Sentence phrase_type;

    typedef TreeGrammar    tree_grammar_type;
    typedef TreeTransducer tree_transducer_type;
    
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    
    typedef HyperGraph     hypergraph_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    
    typedef tree_transducer_type::rule_pair_set_type tree_rule_pair_set_type;
    typedef tree_transducer_type::rule_pair_type     tree_rule_pair_type;
    typedef tree_transducer_type::rule_type          tree_rule_type;
    typedef tree_transducer_type::rule_ptr_type      tree_rule_ptr_type;
    
    typedef std::vector<tree_transducer_type::edge_type, std::allocator<tree_transducer_type::edge_type> > edge_set_type;
    
    typedef std::vector<tree_transducer_type::id_type, std::allocator<tree_transducer_type::id_type> > node_queue_type;
    
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > frontier_type;
    typedef std::deque<frontier_type, std::allocator<frontier_type> > frontier_queue_type;
    
    typedef std::deque<feature_set_type, std::allocator<feature_set_type> >    feature_queue_type;
    typedef std::deque<attribute_set_type, std::allocator<attribute_set_type> > attribute_queue_type;

    typedef feature_set_type::feature_type     feature_type;
    typedef attribute_set_type::attribute_type attribute_type;
    
    // for phrasal matching...
    
#ifdef HAVE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_set<phrase_type, boost::hash<phrase_type>,  std::equal_to<phrase_type>, std::allocator<phrase_type> > phrase_set_type;
#else
    typedef sgi::hash_set<phrase_type, boost::hash<phrase_type>,  std::equal_to<phrase_type>, std::allocator<phrase_type> > phrase_set_type;
#endif
    typedef std::vector<phrase_set_type, std::allocator<phrase_set_type> > phrase_map_type;

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
    
    struct State
    {
      frontier_type            frontier;
      feature_set_type         features;
      attribute_set_type       attributes;
      tree_transducer_type::id_type node;
      
      State() : frontier(), features(), attributes(), node(tree_transducer_type::id_type(-1)) {}
      State(const frontier_type& __frontier,
	    const tree_transducer_type::id_type& __node)
	: frontier(__frontier), features(), attributes(), node(__node) {}
      State(const frontier_type& __frontier,
	    const feature_set_type& __features,
	    const attribute_set_type& __attributes,
	    const tree_transducer_type::id_type& __node)
      : frontier(__frontier), features(__features), attributes(__attributes), node(__node) {}
    };
    typedef State state_type;

    typedef std::deque<state_type, std::allocator<state_type> > queue_type;

    struct NodeMap
    {
      typedef google::dense_hash_map<symbol_type, hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<symbol_type> > node_map_type;

      typedef node_map_type::value_type value_type;
      
      typedef node_map_type::iterator       iterator;
      typedef node_map_type::const_iterator const_iterator;

      NodeMap() : node_map() { node_map.set_empty_key(symbol_type()); }
      
      std::pair<iterator, bool> insert(const value_type& x) { return node_map.insert(x); }
      
      const_iterator find(const symbol_type& x) const { return node_map.find(x); }
      iterator find(const symbol_type& x) { return node_map.find(x); }
      
      const_iterator begin() const { return node_map.begin(); }
      iterator begin() { return node_map.begin(); }

      const_iterator end() const { return node_map.end(); }
      iterator end() { return node_map.end(); }

      node_map_type node_map;
    };

    typedef NodeMap node_map_type;
    typedef std::vector<node_map_type, std::allocator<node_map_type> > node_map_set_type;

    typedef Semiring semiring_type;
    typedef Semiring score_type;
    
    typedef Function function_type;
    
    struct TreeCandidate
    {
      score_type    score;
      
      tree_rule_ptr_type rule;
      
      feature_set_type   features;
      attribute_set_type attributes;
      
      TreeCandidate() : score(), rule(), features(), attributes() {}
      TreeCandidate(const score_type& __score, const tree_rule_ptr_type& __rule, const feature_set_type& __features, const attribute_set_type& __attributes)
	: score(__score), rule(__rule), features(__features), attributes(__attributes) {}
    };
    
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
    
    typedef TreeCandidate tree_candidate_type;
    typedef RuleCandidate rule_candidate_type;
    
    typedef utils::chunk_vector<tree_candidate_type, 4096 / sizeof(tree_candidate_type), std::allocator<tree_candidate_type> > tree_candidate_set_type;
    typedef utils::chunk_vector<rule_candidate_type, 4096 / sizeof(rule_candidate_type), std::allocator<rule_candidate_type> > rule_candidate_set_type;
    

    typedef std::vector<const rule_candidate_type*, std::allocator<const rule_candidate_type*> > rule_candidate_ptr_set_type;
    typedef std::vector<const tree_candidate_type*, std::allocator<const tree_candidate_type*> > tree_candidate_ptr_set_type;
    
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<tree_transducer_type::id_type, tree_candidate_ptr_set_type, utils::hashmurmur<size_t>, std::equal_to<tree_transducer_type::id_type>,
				    std::allocator<std::pair<const tree_transducer_type::id_type, tree_candidate_ptr_set_type> > > tree_candidate_map_type;
    typedef std::tr1::unordered_map<transducer_type::id_type, rule_candidate_ptr_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
				    std::allocator<std::pair<const transducer_type::id_type, rule_candidate_ptr_set_type> > > rule_candidate_map_type;
#else
    typedef sgi::hash_map<tree_transducer_type::id_type, tree_candidate_ptr_set_type, utils::hashmurmur<size_t>, std::equal_to<tree_transducer_type::id_type>,
			  std::allocator<std::pair<const tree_transducer_type::id_type, tree_candidate_ptr_set_type> > > tree_candidate_map_type;
    typedef sgi::hash_map<transducer_type::id_type, rule_candidate_ptr_set_type, utils::hashmurmur<size_t>, std::equal_to<transducer_type::id_type>,
			  std::allocator<std::pair<const transducer_type::id_type, rule_candidate_ptr_set_type> > > rule_candidate_map_type;
#endif
    typedef std::vector<tree_candidate_map_type, std::allocator<tree_candidate_map_type> > tree_candidate_table_type;
    typedef std::vector<rule_candidate_map_type, std::allocator<rule_candidate_map_type> > rule_candidate_table_type;

    struct Candidate
    {      
      score_type    score;

      typename tree_candidate_ptr_set_type::const_iterator tree_first;
      typename tree_candidate_ptr_set_type::const_iterator tree_last;
      
      typename rule_candidate_ptr_set_type::const_iterator rule_first;
      typename rule_candidate_ptr_set_type::const_iterator rule_last;
      
      frontier_type frontier;
      
      feature_set_type   features;
      attribute_set_type attributes;

      bool is_rule() const { return rule_first != rule_last; }
      bool is_tree() const { return rule_first == rule_last; }
      
      score_type candidate_score() const
      {
	return (rule_first == rule_last ? (*tree_first)->score : (*rule_first)->score) * score;
      }
      
      Candidate()
	: score(), tree_first(), tree_last(), rule_first(), rule_last(), frontier(), features(), attributes() { rule_first = rule_last; }
      Candidate(const score_type& __score,
		typename tree_candidate_ptr_set_type::const_iterator __tree_first,
		typename tree_candidate_ptr_set_type::const_iterator __tree_last,
		const frontier_type& __frontier,
		const feature_set_type& __features,
		const attribute_set_type& __attributes)
	: score(__score),
	  tree_first(__tree_first), tree_last(__tree_last), rule_first(), rule_last(),
	  frontier(__frontier),
	  features(__features),
	  attributes(__attributes) { rule_first = rule_last; }
      
      Candidate(typename rule_candidate_ptr_set_type::const_iterator __rule_first,
		typename rule_candidate_ptr_set_type::const_iterator __rule_last)
	: score(semiring::traits<score_type>::one()),
	  tree_first(), tree_last(), rule_first(__rule_first), rule_last(__rule_last),
	  frontier(),
	  features(),
	  attributes() { tree_first = tree_last; }
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

    typedef std::vector<score_type,  std::allocator<score_type> >  score_set_type;
    
    
    ParseTree(const symbol_type& __goal,
	      const tree_grammar_type& __tree_grammar,
	      const grammar_type& __grammar,
	      const function_type& __function,
	      const int __beam_size,
	      const bool __yield_source)
      : goal(__goal),
	tree_grammar(__tree_grammar), 
	grammar(__grammar),
	function(__function),
	beam_size(__beam_size),
	yield_source(__yield_source),
	attr_source_root("source-root")
    {  
      goal_rule = rule_type::create(rule_type(vocab_type::GOAL,
					      rule_type::symbol_set_type(1, goal.non_terminal())));
    }
    
    void operator()(const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      graph_out.clear();
      node_map.clear();
      
      if (! graph_in.is_valid()) return;
      
      node_map.reserve(graph_in.nodes.size());
      node_map.resize(graph_in.nodes.size());
      
      phrase_map.clear();
      phrase_map.reserve(graph_in.nodes.size());
      phrase_map.resize(graph_in.nodes.size());

      label_map.clear();
      
      rule_candidates.clear();
      tree_candidates.clear();
      
      rule_tables.clear();
      rule_tables.reserve(grammar.size());
      rule_tables.resize(grammar.size());
      
      tree_tables.clear();
      tree_tables.reserve(tree_grammar.size());
      tree_tables.resize(tree_grammar.size());
      
      scores.clear();
      scores.reserve(graph_in.nodes.size());
      scores.resize(graph_in.nodes.size(), semiring::traits<score_type>::min());
      
      // bottom-up topological order
      for (size_t id = 0; id != graph_in.nodes.size(); ++ id) {
	candidates.clear();
	heap.clear();
	
	match_tree(id, graph_in);

	// we will do exact phrase-matching...
	if (! grammar.empty())
	  match_phrase(id, graph_in, graph_out);
	
	enumerate_tree(id, graph_in, graph_out);
      }
      
      // goal node... the goal must be mapped goal...
      typename node_map_type::const_iterator niter = node_map[graph_in.goal].find(goal.non_terminal());
      if (niter != node_map[graph_in.goal].end()) {
	// goal node...
	graph_out.goal = graph_out.add_node().id;
	
	// hyperedge to goal...
	hypergraph_type::edge_type& edge = graph_out.add_edge(&(niter->second), (&(niter->second)) + 1);
	edge.rule = goal_rule;
	
	// connect!
	graph_out.connect_edge(edge.id, graph_out.goal);
      }
      
      node_map.clear();
      
      if (graph_out.is_valid())
	graph_out.topologically_sort();
    }

  private:

    void enumerate_tree(const int id, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      rules_enumerated.clear();

      for (int num_pop = 0; ! heap.empty() && num_pop != beam_size; ++ num_pop) {
	const candidate_type* item = heap.top();
	heap.pop();
	
	if (item->is_tree()) {
	  const tree_candidate_type& rule = *(*(item->tree_first));
	  
	  scores[id] = std::max(scores[id], item->score * rule.score);
	  
	  apply_rule(*rule.rule,
		     id,
		     item->frontier,
		     rule.features + item->features,
		     rule.attributes + item->attributes,
		     graph_in,
		     graph_out);
	  
	  // next queue!
	  ++ const_cast<candidate_type*>(item)->tree_first;
	  if (item->tree_first != item->tree_last)
	    heap.push(item);
	} else {
	  const rule_candidate_type& rule = *(*(item->rule_first));
	  
	  scores[id] = std::max(scores[id], item->score * rule.score);

	  rules_enumerated.push_back(*(item->rule_first));
	  
	  ++ const_cast<candidate_type*>(item)->rule_first;
	  if (item->rule_first != item->rule_last)
	    heap.push(item);
	}
      }

      if (! rules_enumerated.empty()) {
	const symbol_type& root_label = graph_in.edges[graph_in.nodes[id].edges.front()].rule->lhs;
	
	typename rule_candidate_ptr_set_type::const_iterator piter_end = rules_enumerated.end();
	for (typename rule_candidate_ptr_set_type::const_iterator piter = rules_enumerated.begin(); piter != piter_end; ++ piter) {
	  const rule_candidate_type& cand = *(*piter);
	  
	  if (node_map[id].find(cand.rule->lhs) == node_map[id].end()) {
	    typename node_map_type::const_iterator liter_end = node_map[id].end();
	    for (typename node_map_type::const_iterator liter = node_map[id].begin(); liter != liter_end; ++ liter) {
	      hypergraph_type::edge_type& edge = graph_out.add_edge();
	      
	      edge.rule = rule_type::create(rule_type(liter->first, cand.rule->rhs));
	      edge.features = cand.features;
	      edge.attributes = cand.attributes;
	      edge.attributes[attr_source_root] = static_cast<const std::string&>(root_label);
	      
	      graph_out.connect_edge(edge.id, liter->second);
	    }
	  }
	  
	  hypergraph_type::edge_type& edge = graph_out.add_edge();
	  
	  edge.rule = cand.rule;
	  edge.features = cand.features;
	  edge.attributes = cand.attributes;
	  edge.attributes[attr_source_root] = static_cast<const std::string&>(root_label);
	  
	  std::pair<typename node_map_type::iterator, bool> result = node_map[id].insert(std::make_pair(edge.rule->lhs, 0));
	  if (result.second)
	    result.first->second = graph_out.add_node().id;
	  
	  graph_out.connect_edge(edge.id, result.first->second);
	}
      }
    }
    
    void match_phrase(const int id, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      typedef std::deque<phrase_type, std::allocator<phrase_type> >  buffer_type;
      
      if (graph_in.nodes[id].edges.empty()) return;
      
      // first, construct prases

      buffer_type buffer;
      buffer_type buffer_next;
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = graph_in.nodes[id].edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = graph_in.nodes[id].edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = graph_in.edges[*eiter];
	
	buffer.clear();
	buffer.push_back(phrase_type());

	int non_terminal_pos = 0;
	rule_type::symbol_set_type::const_iterator titer_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator titer = edge.rule->rhs.begin(); titer != titer_end; ++ titer)
	  if (titer->is_non_terminal()) {
	    const int __non_terminal_index = titer->non_terminal_index();
	    const int pos = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	    ++ non_terminal_pos;
	    
	    // combine buffer and tails...
	    buffer_next.clear();
	    
	    phrase_set_type::const_iterator piter_end = phrase_map[edge.tails[pos]].end();
	    for (phrase_set_type::const_iterator piter = phrase_map[edge.tails[pos]].begin(); piter != piter_end; ++ piter) {
	      buffer_type::const_iterator biter_end = buffer.end();
	      for (buffer_type::const_iterator biter = buffer.begin(); biter != biter_end; ++ biter) {
		buffer_next.push_back(*biter);
		buffer_next.back().insert(buffer_next.back().end(), piter->begin(), piter->end());
	      }
	    }
	    
	    buffer.swap(buffer_next);
	  } else if (*titer != vocab_type::EPSILON) {
	    buffer_type::iterator biter_end = buffer.end();
	    for (buffer_type::iterator biter = buffer.begin(); biter != biter_end; ++ biter)
	      biter->push_back(*titer);
	  }
	
	phrase_map[id].insert(buffer.begin(), buffer.end());
      }
      
      // then, try matching within this span...

      for (size_t grammar_id = 0; grammar_id != grammar.size(); ++ grammar_id) {
	const transducer_type& transducer = grammar[grammar_id];
	
	phrase_set_type::const_iterator piter_end = phrase_map[id].end();
	for (phrase_set_type::const_iterator piter = phrase_map[id].begin(); piter != piter_end; ++ piter) {
	  const phrase_type& phrase = *piter;
	  
	  transducer_type::id_type node = transducer.root();
	  
	  phrase_type::const_iterator iter_end = phrase.end();
	  for (phrase_type::const_iterator iter = phrase.begin(); iter != iter_end; ++ iter) {
	    node = transducer.next(node, *iter);
	    if (node == transducer.root()) break;
	  }
	  
	  if (node == transducer.root()) continue;
	  
	  const rule_candidate_ptr_set_type& rules = candidate_rules(grammar_id, node);
	  
	  if (rules.empty()) continue;
	  
	  candidates.push_back(candidate_type(rules.begin(), rules.end()));
	  heap.push(&candidates.back());
	}
      }
    }
    
    void match_tree(const int id, const hypergraph_type& graph_in)
    {
      if (graph_in.nodes[id].edges.empty()) return;
      
      //std::cerr << "node: " << id << std::endl;
      
      queue_type queue;
      
      for (size_t grammar_id = 0; grammar_id != tree_grammar.size(); ++ grammar_id) {
	const tree_transducer_type& transducer = tree_grammar[grammar_id];

	//std::cerr << "transducer: " << grammar_id << std::endl;
	
	const symbol_type cat = graph_in.edges[graph_in.nodes[id].edges.front()].rule->lhs;
	const tree_transducer_type::edge_type edge_id = transducer.edge(cat);
	
	if (edge_id == tree_transducer_type::edge_type()) continue;
	
	const tree_transducer_type::edge_type edge_epsilon = transducer.edge(vocab_type::EPSILON);
	const tree_transducer_type::edge_type edge_none    = transducer.edge(vocab_type::NONE);
	
	tree_transducer_type::id_type node = transducer.next(transducer.root(), edge_id);
	if (node == transducer.root()) continue;
	
	node = transducer.next(node, edge_none);
	if (node == transducer.root()) continue;
	
	//std::cerr << "grammar cat: " << cat << " id: " << edge_id << std::endl;
	
	queue.clear();
	queue.push_back(state_type(frontier_type(1, id), node));
	
	while (! queue.empty()) {
	  const state_type& state = queue.front();
	  
	  frontier_queue_type frontiers(1, frontier_type());
	  frontier_queue_type frontiers_next;
	  
	  node_queue_type nodes(1, state.node);
	  node_queue_type nodes_next;

	  feature_queue_type features(1, state.features);
	  feature_queue_type features_next;
	  
	  attribute_queue_type attributes(1, state.attributes);
	  attribute_queue_type attributes_next;

	  edge_set_type edges;
	  
	  frontier_type::const_iterator niter_end = state.frontier.end();
	  for (frontier_type::const_iterator niter = state.frontier.begin(); niter != niter_end; ++ niter) {
	    
	    frontiers_next.clear();
	    nodes_next.clear();
	    features_next.clear();
	    attributes_next.clear();
	    
	    edges.clear();
	    hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = graph_in.nodes[*niter].edges.end();
	    for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = graph_in.nodes[*niter].edges.begin(); eiter != eiter_end; ++ eiter) {
	      const hypergraph_type::edge_type& edge = graph_in.edges[*eiter];
	      
	      edges.push_back(transducer.edge(edge.rule->rhs));
	    }
	    
	    frontier_queue_type::const_iterator  fiter = frontiers.begin();
	    feature_queue_type::const_iterator   siter = features.begin();
	    attribute_queue_type::const_iterator aiter = attributes.begin();

	    node_queue_type::const_iterator titer_end = nodes.end();
	    for (node_queue_type::const_iterator titer = nodes.begin(); titer != titer_end; ++ titer, ++ fiter, ++ siter, ++ aiter) {
	      const tree_transducer_type::id_type node_epsilon = transducer.next(*titer, edge_epsilon);
	      if (node_epsilon != transducer.root()) {
		frontier_type frontier(*fiter);
		frontier.push_back(*niter);
		
		frontiers_next.push_back(frontier);
		nodes_next.push_back(node_epsilon);
		features_next.push_back(*siter);
		attributes_next.push_back(*aiter);
	      }

	      edge_set_type::const_iterator eiter_begin = edges.begin();
	      edge_set_type::const_iterator eiter_end = edges.end();
	      for (edge_set_type::const_iterator eiter = eiter_begin; eiter != eiter_end; ++ eiter)
		if (*eiter != tree_transducer_type::edge_type()) {
		  const tree_transducer_type::edge_type& edge_id = *eiter;
		  
		  const tree_transducer_type::id_type node_edge = transducer.next(*titer, edge_id);
		  if (node_edge != transducer.root()) {
		    const hypergraph_type::edge_type& edge = graph_in.edges[graph_in.nodes[*niter].edges[eiter - eiter_begin]];
		    
		    frontier_type frontier(*fiter);
		    frontier.insert(frontier.end(), edge.tails.begin(), edge.tails.end());
		    
		    frontiers_next.push_back(frontier);
		    nodes_next.push_back(node_edge);
		    features_next.push_back(*siter + edge.features);
		    attributes_next.push_back(*aiter + edge.attributes);
		  }
		}
	    }
	    
	    frontiers.swap(frontiers_next);
	    nodes.swap(nodes_next);
	    features.swap(features_next);
	    attributes.swap(attributes_next);
	    
	    frontiers_next.clear();
	    nodes_next.clear();
	    features_next.clear();
	    attributes_next.clear();
	  }

	  //std::cerr << "finished loop: " << frontiers.size() << std::endl;
	  
	  // frontiers and nodes contain new frontier!
	  // in addition, we need to traverse transducer.next() with edge_none!
	  
	  frontier_queue_type::const_iterator  fiter = frontiers.begin();
	  feature_queue_type::const_iterator   siter = features.begin();
	  attribute_queue_type::const_iterator aiter = attributes.begin();
	  
	  node_queue_type::const_iterator titer_end = nodes.end();
	  for (node_queue_type::const_iterator titer = nodes.begin(); titer != titer_end; ++ titer, ++ fiter, ++ siter, ++ aiter) {
	    const tree_transducer_type::id_type node_none = transducer.next(*titer, edge_none);
	    if (node_none == transducer.root()) continue;
	    
	    queue.push_back(state_type(*fiter, node_none));
	    
	    const tree_candidate_ptr_set_type& rules = candidate_trees(grammar_id, node_none);
	    
	    if (rules.empty()) continue;
	    
	    // compute frontier scores
	    score_type score = semiring::traits<score_type>::one();
	    frontier_type::const_iterator iter_end = fiter->end();
	    for (frontier_type::const_iterator iter = fiter->begin(); iter != iter_end; ++ iter)
	      score *= scores[*iter];
	    
	    candidates.push_back(candidate_type(score, rules.begin(), rules.end(), *fiter, *siter, *aiter));
	    heap.push(&candidates.back());
	  }
	  
	  queue.pop_front();
	}
      }
    }
    
    
    void apply_rule(const tree_rule_type& rule,
		    const hypergraph_type::id_type root_in,
		    const frontier_type& frontiers,
		    const feature_set_type& features,
		    const attribute_set_type& attributes,
		    const hypergraph_type& graph_in,
		    hypergraph_type& graph_out)
    {
      // apply rule for the frontiers with root at graph_in and generate graph_out...
      // frontiers are ordered wrt source side rule's frontier nodes (in pre-traversal order),## the same as post-traversal order...
      //
      // collect frontiers...
      
      //
      // construct graph_out in pre-order...
      //

      const symbol_type& root_label = graph_in.edges[graph_in.nodes[root_in].edges.front()].rule->lhs;
      
      std::pair<typename node_map_type::iterator, bool> result = node_map[root_in].insert(std::make_pair(rule.label.non_terminal(), 0));
      if (result.second)
	result.first->second = graph_out.add_node().id;
      
      int non_terminal_pos = 0;
      
      const hypergraph_type::id_type edge_id = construct_graph(rule, result.first->second, frontiers, graph_out, non_terminal_pos);
      
      graph_out.edges[edge_id].features   = features;
      graph_out.edges[edge_id].attributes = attributes;
      
      // root-label is assigned to source-root attribute
      graph_out.edges[edge_id].attributes[attr_source_root] = static_cast<const std::string&>(root_label);
    }
    
    hypergraph_type::id_type construct_graph(const tree_rule_type& rule,
					     hypergraph_type::id_type root,
					     const frontier_type& frontiers,
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
	    
	    const hypergraph_type::id_type node = frontiers[non_terminal_index];
	    
	    std::pair<typename node_map_type::iterator, bool> result = node_map[node].insert(std::make_pair(aiter->label.non_terminal(), 0));
	    if (result.second)
	      result.first->second = graph.add_node().id;
	    
	    tails.push_back(result.first->second);
	  } else
	    tails.push_back(0);
	  
	  rhs.push_back(aiter->label.non_terminal());
	} else
	  rhs.push_back(aiter->label);

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
      } else {
	edge_id = graph.add_edge(tails.begin(), tails.end()).id;
	graph.edges[edge_id].rule = rule_type::create(rule_type(rule.label, rhs.begin(), rhs.end()));
	graph.connect_edge(edge_id, root);
      }
      
      return edge_id;
    }

    template <typename Tp>
    struct greater_ptr_score
    {
      bool operator()(const Tp* x, const Tp* y) const
      {
	return x->score > y->score;
      }
    };

    const rule_candidate_ptr_set_type& candidate_rules(const size_type& table, const transducer_type::id_type& node)
    {
      typename rule_candidate_map_type::iterator riter = rule_tables[table].find(node);
      if (riter == rule_tables[table].end()) {
	riter = rule_tables[table].insert(std::make_pair(node, rule_candidate_ptr_set_type())).first;
	
	const transducer_type::rule_pair_set_type& rules = grammar[table].rules(node);
	
	if (rules.size() > beam_size) {
	  transducer_type::rule_pair_set_type::const_iterator iter_begin = rules.begin();
	  transducer_type::rule_pair_set_type::const_iterator iter_end   = rules.end();
	  for (transducer_type::rule_pair_set_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter) {
	    const score_type score = function(iter->features);
	    
	    if (riter->second.size() < beam_size || score >= riter->second.front()->score) {
	      rule_candidates.push_back(rule_candidate_type(score,
							    yield_source ? iter->source : iter->target,
							    iter->features,
							    iter->attributes));
	      
	      riter->second.push_back(&(rule_candidates.back()));
	      
	      std::push_heap(riter->second.begin(), riter->second.end(), greater_ptr_score<rule_candidate_type>());
	    }
	  }
	  
	  // heap-sort!
	  std::sort_heap(riter->second.begin(), riter->second.end(), greater_ptr_score<rule_candidate_type>());
	  
	  // shrink...
	  rule_candidate_ptr_set_type(riter->second).swap(riter->second);
	} else {
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
      }
      return riter->second;
    }

    const tree_candidate_ptr_set_type& candidate_trees(const size_type& table, const tree_transducer_type::id_type& node)
    {
      typename tree_candidate_map_type::iterator riter = tree_tables[table].find(node);
      if (riter == tree_tables[table].end()) {
	riter = tree_tables[table].insert(std::make_pair(node, tree_candidate_ptr_set_type())).first;
	
	const tree_transducer_type::rule_pair_set_type& rules = tree_grammar[table].rules(node);
	
	if (rules.size() > beam_size) {
	  tree_transducer_type::rule_pair_set_type::const_iterator iter_begin = rules.begin();
	  tree_transducer_type::rule_pair_set_type::const_iterator iter_end   = rules.end();
	  for (tree_transducer_type::rule_pair_set_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter) {
	    const score_type score = function(iter->features);

	    if (riter->second.size() < beam_size || score >= riter->second.front()->score) {
	      tree_candidates.push_back(tree_candidate_type(score,
							    yield_source ? iter->source : iter->target,
							    iter->features,
							    iter->attributes));
	      
	      riter->second.push_back(&(tree_candidates.back()));
	      
	      std::push_heap(riter->second.begin(), riter->second.end(), greater_ptr_score<tree_candidate_type>());
	    }
	  }
	  
	  // heap-sort!
	  std::sort_heap(riter->second.begin(), riter->second.end(), greater_ptr_score<tree_candidate_type>());
	  
	  // shrink...
	  tree_candidate_ptr_set_type(riter->second).swap(riter->second);
	} else {
	  riter->second.reserve(rules.size());
	
	  tree_transducer_type::rule_pair_set_type::const_iterator iter_begin = rules.begin();
	  tree_transducer_type::rule_pair_set_type::const_iterator iter_end   = rules.end();
	  for (tree_transducer_type::rule_pair_set_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter) {
	    tree_candidates.push_back(tree_candidate_type(function(iter->features),
							  yield_source ? iter->source : iter->target,
							  iter->features,
							  iter->attributes));
	  
	    riter->second.push_back(&(tree_candidates.back()));
	  }
	
	  std::sort(riter->second.begin(), riter->second.end(), greater_ptr_score<tree_candidate_type>());
	}
      }
      return riter->second;
    }
    
  private:

    node_map_set_type node_map;
    
    phrase_map_type phrase_map;

    internal_label_map_type label_map;
    
    candidate_set_type    candidates;
    candidate_heap_type   heap;

    rule_candidate_ptr_set_type rules_enumerated;

    rule_candidate_set_type   rule_candidates;
    rule_candidate_table_type rule_tables;
    
    tree_candidate_set_type   tree_candidates;
    tree_candidate_table_type tree_tables;
    
    score_set_type scores;

    rule_ptr_type goal_rule;
    
    symbol_type goal;
    const tree_grammar_type& tree_grammar;
    const grammar_type& grammar;
    const function_type& function;
    
    const int beam_size;
    const bool yield_source;
    
    attribute_type attr_source_root;
  };
  
  template <typename Function>
  inline
  void parse_tree(const Symbol& goal, const TreeGrammar& tree_grammar, const Grammar& grammar, const Function& function, const HyperGraph& graph_in, HyperGraph& graph_out, const int size, const bool yield_source=false)
  {
    ParseTree<typename Function::value_type, Function> __parser(goal, tree_grammar, grammar, function, size, yield_source);
    __parser(graph_in, graph_out);
  }
};

#endif
