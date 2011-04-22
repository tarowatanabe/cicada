// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_TREE__HPP__
#define __CICADA__PARSE_TREE__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>

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
#include <utils/b_heap.hpp>
#include <utils/std_heap.hpp>
#include <utils/bithack.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

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
    
    struct RuleCandidate
    {
      score_type score;
      
      tree_rule_ptr_type rule;
      
      feature_set_type   features;
      attribute_set_type attributes;

      RuleCandidate() : score(), rule(), features(), attributes() {}
      RuleCandidate(const score_type& __score, const tree_rule_ptr_type& __rule, const feature_set_type& __features, const attribute_set_type& __attributes)
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

    struct Candidate
    {
      typename rule_candidate_ptr_set_type::const_iterator first;
      typename rule_candidate_ptr_set_type::const_iterator last;
      
      score_type    score;
      
      frontier_type frontier;
      
      feature_set_type   features;
      attribute_set_type attributes;
      
      score_type candidate_score() const
      {
	return first->score * score;
      }
      
      Candidate() : first(), last(), score(), frontier(), features(), attributes() {}
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
      
      rule_candidates.clear();
      rule_tables.clear();
      rule_tables.reserve(tree_grammar.size());
      rule_tables.resize(tree_grammar.size());
      
      scores_max.clear();
      scores_min.clear();
      scores_max.reserve(graph_in.nodes.size());
      scores_min.reserve(graph_in.nodes.size());
      scores_max.resize(graph_in.nodes.size(), semiring::traits<score_type>::min());
      scores_min.resize(graph_in.nodes.size(), semiring::traits<score_type>::max());
      
      // bottom-up topological order
      for (size_t id = 0; id != graph_in.nodes.size(); ++ id) {
	candidates.clear();
	heap.clear();
	
	match_tree(id, graph_in);

	enumerate_tree(id, graph_in, graph_out);
	
	// we will do exact phrase-matching...
	if (! grammar.empty())
	  match_phrase(id, graph_in, graph_out);
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
    }

  private:

    void enumerate_tree(const int id, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      for (int num_pop = 0; ! heap.empty() && num_pop != beam_size; ++ num_pop) {
	const candidate_type* item = heap.top();
	heap.pop();
	
	const rule_candidate_type& rule = *(*(item->first));
	const score_type score = item->score * rule.score;
	
	// update scores...
	scores_max[id] = std::max(scores_max[id], score);
	scores_min[id] = std::min(scores_min[id], score);
	
	apply_rule(rule.rule,
		   id,
		   item->frontier,
		   rule.features + item->features,
		   rule.attributes + item->attributes,
		   graph_in, graph_out);
	
	// next queue!
	++ const_cast<candidate_type*>(item)->first;
	if (item->first != item->last)
	  heap.push(item);
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
      
      const symbol_type& root_label = graph_in.edges[graph_in.nodes[id].edges.front()].rule->lhs;

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
	  
	  const transducer_type::rule_pair_set_type& rules = transducer.rules(node);
	  
	  if (rules.empty()) continue;
	  
	  transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	  for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	    
	    const score_type score = function(riter->features);
	    
	    // we do not include bad phrases
	    if (score < scores_min[id]) continue;
	    
	    scores_max[id] = std::max(scores_max[id], score);
	    
	    const rule_ptr_type rule = (yield_source ? riter->source : riter->target);
	    
	    if (node_map[id].find(rule->lhs) == node_map[id].end()) {
	      // we will try all the combination of lhs in node_map[id]
	      
	      typename node_map_type::const_iterator liter_end = node_map[id].end();
	      for (typename node_map_type::const_iterator liter = node_map[id].begin(); liter != liter_end; ++ liter) {
		hypergraph_type::edge_type& edge = graph_out.add_edge();
		edge.rule = rule_type::create(rule_type(liter->first, rule->rhs.begin(), rule->rhs.end()));
		edge.features = riter->features;
		edge.attributes = riter->attributes;
		
		edge.attributes[attr_source_root] = static_cast<const std::string&>(root_label);
		
		graph_out.connect_edge(edge.id, liter->second);
	      }
	    }
	    
	    hypergraph_type::edge_type& edge = graph_out.add_edge();
	    edge.rule = rule;
	    edge.features = riter->features;
	    edge.attributes = riter->attributes;
	    
	    edge.attributes[attr_source_root] = static_cast<const std::string&>(root_label);
	    
	    std::pair<typename node_map_type::iterator, bool> result = node_map[id].insert(std::make_pair(edge.rule->lhs, 0));
	    if (result.second)
	      result.first->second = graph_out.add_node().id;
	    
	    graph_out.connect_edge(edge.id, result.first->second);
	  }
	}
      }
    }
    
    void match_tree(const int id, const hypergraph_type& graph_in, hypergraph_type& graph_out)
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
	    
	    const rule_candidate_ptr_set_type& rules = cands(grammar_id, node_none);
	    
	    if (rules.empty()) continue;
	    
	    // compute frontier scores
	    score_type score = semiring::traits<score_type>::one();
	    frontier_type::const_iterator iter_end = fiter->end();
	    for (frontier_type::const_iterator iter = fiter->begin(); iter != iter_end; ++ iter)
	      score *= scores_max[*iter];
	    
	    candidates.push_back(candidate_type());
	    candidate_type& cand = candidates.back();
	    cand.first  = rules.begin();
	    cand.second = rules.end();
	    
	    cand.score = score;
	    
	    cand.frontier   = *fiter;
	    cand.features   = *siter;
	    cand.attributes = *aiter;
	    
	    heap.push(&cand);
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
      
      const hypergraph_type::id_type edge_id = construct_graph(rule, result.first->second, frontiers, graph_in, graph_out, non_terminal_pos);
      
      graph_out.edges[edge_id].features   = features;
      graph_out.edges[edge_id].attributes = attributes;
      
      // root-label is assigned to source-root attribute
      graph_out.edges[edge_id].attributes[attr_source_root] = static_cast<const std::string&>(root_label);
    }
    
    hypergraph_type::id_type construct_graph(const tree_rule_type& rule,
					     const hypergraph_type::id_type root,
					     const frontier_type& frontiers,
					     const hypergraph_type& graph_in,
					     hypergraph_type& graph_out,
					     int& non_terminal_pos)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > rhs_type;
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails_type;
      
      rhs_type rhs;
      tails_type tails;
      
      int pos = 1;
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
	      result.first->second = graph_out.add_node().id;
	    
	    tails.push_back(result.first->second);
	  } else
	    tails.push_back(graph_out.add_node().id);
	  
	  rhs.push_back(aiter->label.non_terminal(pos));
	  
	  ++ pos;
	} else
	  rhs.push_back(aiter->label);
      
      const hypergraph_type::id_type edge_id = graph_out.add_edge(tails.begin(), tails.end()).id;
      graph_out.edges[edge_id].rule = rule_type::create(rule_type(rule.label, rhs.begin(), rhs.end()));
      
      //std::cerr << "output forest rule: " << *(graph_out.edges[edge_id].rule) << std::endl;
      
      graph_out.connect_edge(edge_id, root);
      
      tails_type::const_iterator titer = tails.begin();
      for (tree_rule_type::const_iterator aiter = rule.begin(); aiter != aiter_end; ++ aiter) {
	if (aiter->label.is_non_terminal() && ! aiter->antecedents.empty()) {
	  construct_graph(*aiter, *titer, frontiers, graph_in, graph_out, non_terminal_pos);
	  ++ titer;
	}
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
    
    const rule_candidate_ptr_set_type& cands(const int grammar_id, const tree_transducer_type::id_type& node)
    {
      typename rule_candidate_map_type::iterator riter = rule_tables[grammar_id].find(node);
      if (riter == rule_tables[grammar_id].end()) {
	riter = rule_tables[grammar_id].insert(std::make_pair(node, rule_candidate_ptr_set_type())).first;
	
	const tree_transducer_type::rule_pair_set_type& rules = tree_grammar[grammar_id].rules(node);

	riter->second.reserve(rules.size());
	
	tree_transducer_type::rule_pair_set_type::const_iterator iter_begin = rules.begin();
	tree_transducer_type::rule_pair_set_type::const_iterator iter_end   = rules.end();
	for (tree_transducer_type::rule_pair_set_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter) {
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

    node_map_set_type node_map;
    
    phrase_map_type phrase_map;
    
    candidate_set_type    candidates;
    candidate_heap_type   heap;
    
    rule_candidate_set_type   rule_candidates;
    rule_candidate_table_type rule_tables;
    
    score_set_type        scores_max;
    score_set_type        scores_min;

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
