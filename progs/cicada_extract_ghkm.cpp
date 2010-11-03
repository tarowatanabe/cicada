
#include <set>

#include "cidada/symbol.hpp"
#include "cicada/alignment.hpp"
#inlcude "cicada/hypergraph.hpp"
#include "cicada/sentence.hpp"

typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Symbol     word_type;
typedef cicada::Symbol     symbol_type;
typedef cicada::Vocab      vocab_type;
typedef cicada::Sentence   sentence_type;
typedef cicada::Alignment  alignment_type;
typedef cicada::TreeRule   tree_rule_type;
typedef cicada::Rule       rule_type;

typedef alignment_type::point_type point_type;

struct RulePair
{
  tree_rule_type source;
  tree_rule_type target;
  alignment_type alignment;
  double count;
  
  RulePair() : source(), target(), alignment(), count(0) {}
  RulePair(const tree_rule_type& __source, const tree_rule_type& __target, const alignment_type& __alignment, const double& __count)
    : source(__source), target(__target), alignment(__alignment), count(__count) {}
};

typedef RulePair rule_pair_type;
typedef std::vector<rule_pair_type, std::allocator<rule_pair_type> > rule_pair_set_type;


//
// forest-based translation rule extraction based on
//
// @InProceedings{mi-huang:2008:EMNLP,
//   author    = {Mi, Haitao  and  Huang, Liang},
//   title     = {Forest-based Translation Rule Extraction},
//   booktitle = {Proceedings of the 2008 Conference on Empirical Methods in Natural Language Processing},
//   month     = {October},
//   year      = {2008},
//   address   = {Honolulu, Hawaii},
//   publisher = {Association for Computational Linguistics},
//   pages     = {206--214},
//   url       = {http://www.aclweb.org/anthology/D08-1022}
// }


// we assume that the hypergraph is parse forest
struct ExtractGHKM
{

  struct Span
  {
    typedef std::set<int, std::less<int>, std::allocator<int> > span_type;

    typedef span_type::const_iterator const_iterator;
    typedef span_type::iterator       iterator;
  
    Span() : span() {}
    Span(const int first, const int last) : span()
    {
      for (int i = first; i < last; ++ i)
	span.insert(i);
    }

    Span& operator|=(const Span& x)
    {
      span.insert(x.span.begin(), x.span.end());
    }

    const_iterator begin() const { return span.begin(); }
    const_iterator end() const { return span.end(); }

    const_iterator lower_bound(int pos) const { return span.lower_bound(pos); }
    
    std::pair<int, int> range() const
    {
      if (span.empty())
	return std::make_pair(0, 0);
      else
	return std::make_pair(*(span.begin()), *(-- span.end()) + 1);
    }
  
    void set(int i)
    {
      span.insert(i);
    }
    
    void clear() { span.clear(); }

    bool intersect(const std::pair<int, int>& range) const
    {
      if (range.first == range.last)
	return false;
      
      if (span.empty())
	return false;
      else
	return span.lower_bound(range.first) != span.lower_bound(range.last);
    }

    bool empty() const { return span.empty(); }
  
    span_type span;
  };

  typedef Span span_type;
  typedef std::vector<span_type, std::allocator<span_type> > span_set_type;

  typedef std::vector<bool, std::allocator<bool> > admissible_set_type;
  
  void operator()(const hypergraph_type& graph,
		  const sentence_type& sentence,
		  const alignment_type& alignment,
		  rule_pair_set_type& rules) const
  {
    spans.clear();
    complements.clear();
    unaligned.clear();
    admissibles.clear();
    derivations.clear();
    node_map.clear();
    
    spans.resize(graph.nodes.size());
    complements.resize(graph.nodes.size());
    unaligned.resize(graph.nodes.size());
    node_map.resize(graph.nodes.size());
    
    admissible_nodes(graph, sentence, alignment);
    
    minimal_rules(graph);
    
  }
  
  void minimal_rules(const hypergraph_type& graph,
		     const sentence_type& sentence)
  {
    // construc derivations wrt non-aligned words...

    size_t goal_node = size_t(-1);
    
    for (size_t id = 0; id != graph.nodes.size(); ++ id) 
      if (admissibles[id]) {
	const hypergraph_type::node_type& node = graph.nodes[id];

	//
	// compute minimal range and maximum outer-range
	//
	
	const range_type range_min(spans[id].range());
	
	span_type::const_iterator riter_first = complements[id].lower_bound(range_min.first);
	span_type::const_iterator riter_last  = complements[id].lower_bound(range_min.second);
	
	const range_type range_max(riter_first == complements[id].begin() ? 0 : *(-- riter_first) + 1,
				   riter_last != complements[id].end() ? *riter_last : static_cast<int>(sentence.size()));
	
	const bool is_goal(id == graph.goal);
	
	range_node_map_type buf;
	buf.set_empty_key(range_type(0, 0));
	
	queue_type queue;
	
	// construct initial frontiers...
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  queue.resize(queue.size() + 1);
	  frontier_type& frontier = queue.back();
	  
	  frontier.first.push_back(*eiter);
	  
	  hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	  for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	    if (! admissible[*titer])
	      frontier.second.push_back(*titer);
	}

	while (! queue.empty()) {
	  const frontier_type& frontier = queue.front();
	  
	  if (frontier.second.empty()) {
	    // construct rule from "fragment", frontier.first
	    // by enumerating edge and its tails in depth-first manner...
	    // use recursive call for simplicity...
	    
	    // compute "tails" for the fragment
	    
	    node_set_type tails;
	    edge_set_type::const_iterator eiter_end = frontier.first.end();
	    for (edge_set_type::const_iterator eiter = frontier.first.begin(); eiter != eiter_end; ++ eiter) {
	      const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	      
	      hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	      for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
		if (admissibles[*titer])
		  tails.push_back(*titer);
	    }
	    
	    // iterate over tails and compute span...
	    // code taken from apply_exact...
	    
	    index_set_type j_ends(tails.size(), 0);
	    index_set_type j(tails.size(), 0);
	    
	    for (size_t i = 0; i != tails.size(); ++ i)
	      j_ends[i] = node_map[tails[i]].size();
	    
	    node_set_type tails_next(tails.size());
	    
	    for (;;) {
	      bool is_valid = true;
	      
	      range_type range(range_min);
	      range_set_type ranges;
	      for (size_t i = 0; is_valid && i != edge.tails.size(); ++ i) {
		tails_next[i] = node_map[tails[i]][j[i]];
		
		const range_type& range_antecedent = derivations[tails_next[i]].range;
		
		if (! ranges.empty()) {
		  range_set_type::const_iterator riter_end = ranges.end();
		  for (range_set_type::const_iterator riter = ranges.begin(); is_valid && riter != riter_end; ++ riter)
		    is_valid = range_antecedent.second <= riter->first || riter->second <= range_antecedent.first;
		}
		
		if (is_valid) {
		  ranges.insert(range_antecedent);
		  
		  range.first  = utils::bithack::min(range.first, range_antecedent.first);
		  range.second = utils::bithack::max(range.second, range_antecedent.second);
		}
	      }
	      
	      if (is_valid) {
		// our tail-ranges are valid... and we will consider alternatives from range to range_max...
		
		if (is_goal) {
		  if (goal_node == size_t(-1)) {
		    goal_node = derivations.size();
		    derivations.resize(goal_node + 1);
		    
		    // range is always [0, sentence.size())
		    derivations.back().range = range_type(0, sentence.size());
		  }
		  
		  derivations[goal_node].edges.push_back(derivation_type(frontier.first, tails_next));
		} else {
		  // we will compute all possible ranges...
		  for (int first = range.first; first <= range_max.first; ++ first)
		    for (int last = range.second; last <= range_max.second; ++ last) {
		      const range_type range_next(first, last);
		      
		      range_node_map_type::iterator biter = buf.find(range_next);
		      if (biter == buf.end()) {
			derivatoins.resize(derivations.size() + 1);
			derivations.back().range = range_next;
			
			node_map[id].push_back(derivations.size() - 1);
			
			biter = buf.insert(std::make_pair(range_next, derivations.size() - 1)).first;
		      }
		      
		      derivations[biter->second].edges.push_back(derivation_type(frontier.first, tails_next));
		    }
		}
	      }
	      
	      size_t index = 0;
	      for (/**/; index != tails.size(); ++ index) {
		++ j[index];
		if (j[index] < j_ends[index]) break;
		j[index] = 0;
	      }
	      
	      // finished!
	      if (index == tails.size()) break;
	    }
	    
	  } else {
	    // incomplete... futher expand!
	    const hypergraph_type::node_type& node = graph.nodes[frontier.second.front()];
	    
	    hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	      const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	      
	      queue.resize(queue.size() + 1);
	      frontier_type& frontier_next = queue.back();
	      
	      frontier_next.first = frontier.first;
	      frontier_next.first.push_back(*eiter);
	      
	      hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	      for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
		if (! admissible[*titer])
		  frontier_next.second.push_back(*titer);
	      
	      frontier_next.second.insert(frontier_next.second.end(), frontier.second.begin() + 1, frontier.second.end());
	    }
	  }
	  
	  queue.pop_front();
	}
      }
  }
  
  void admissible_nodes(const hypergraph_type& graph,
			const sentence_type& sentence,
			const alignment_type& alignment)
  {
    typedef std::multimap<int, int, std::less<int>, std::allocator<std::pair<const int, int> > > alignment_map_type;
    
    alignment_map_type aligns(alignment.begin(), alignment.end());
    
    // first, bottom-up traversal to compute span...
    for (size_t id = 0; id != graph.nodes.size(); ++ id) {
      const hypergraph_type::node_type& node = graph.nodes[id];
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	
	alignment_map_type::const_iterator liter = alignment_map.lower_bound(edge.first);
	alignment_map_type::const_iterator uiter = alignment_map.lower_bound(edge.last);
	
	for (alignment_map_type::const_iterator iter = liter; iter != uiter; ++ iter)
	  spans[id].set(iter->target);
      }
    }
    
    // second, top-down traversal to compute complement
    for (int id = graph.nodes.size(); id >= 0; -- id) {
      const hypergraph_type::node_type& node = graph.nodes[id];
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	
	hypergraph_type::edge_type::node_set_type::const_iterator titer_begin = edge.tails.begin();
	hypergraph_type::edge_type::node_set_type::const_iterator titer_end   = edge.tails.end();
	
	for (hypergraph_type::edge_type::node_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	  complements[*titer] |= complements[id];
	  for (hypergraph_type::edge_type::node_set_type::const_iterator iiter = titer_begin; iiter != titer_end; ++ iiter)
	    if (iiter != titer)
	      complements[*titer] |= spans[*iiter];
	}
      }
    }

    // check whether admissible or not...
    // do we also compute non-aligned words...?
    for (size_t id = 0; id != graph.nodes.size(); ++ id)
      admissibles[id] = ! spans[id].empty() && ! complements[id].intersect(spans[id].range());
  }
  
  span_set_type spans;
  span_set_type complements;
  span_set_type unaligned;

  admissible_set_type admissibles;

  derivation_graph_type derivations;
};
