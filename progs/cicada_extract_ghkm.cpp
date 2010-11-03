
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
    
    spans.resize(graph.nodes.size());
    complements.resize(graph.nodes.size());
    unaligned.resize(graph.nodes.size());
    admissibles.resize(graph.nodes.size());
    
    admissible_nodes(graph, sentence, alignment);
    
    
    minimal_rules(graph);
  }
  
  void minimal_rules(const hypergraph_type& graph)
  {
    // How do we attach non-aligned words...?
    
    for (sizse_t id = 0; id != graph.nodes.size(); ++ id) 
      if (admissibles[id]) {
	const hypergraph_type::node_type& node = graph.nodes[id];
	
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
	    
	    // construct derivation graph... ignoring non-aligned words...
	    // how to handle non-aligned words...?
	    
	    
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
    for (size_t id = 0; id != graph.nodes.size(); ++ id) {
      const std::pair<int, int> range = spans[id].range();
      
      admissibles[id] = ! spans[id].empty() && ! complements[id].intersect(range);
      
      if (! admissibles[id]) continue;
      
      // compute unaliged....
      
      span_type::const_iterator citer_first = complements[id].lower_bound(range.first);
      if (citer_first != complements[id].begin())
	-- citer_first;
      const int first = (citer_first != complements[id].end() ? *citer_first + 1: 0);
      for (int i = first; i < range.first; ++ i)
	unaligned[id].insert(i);
      
      span_type::const_iterator siter_prev = spans[id].begin();
      span_type::const_iterator siter_end = spans[id].end();
      span_type::const_iterator siter = spans[id].begin();
      ++ siter;
      
      for (/**/; siter != siter_end; ++ siter) {
	for (int i = *siter_prev + 1; i != *siter; ++ i)
	  unaligned[id].insert(i);
	siter_prev = siter;
      }
      
      span_type::const_iterator citer_last = complements[id].lower_bound(range.second);
      const int last = (citer_last != complements[id].end() ? *citer_last : sentence.size());
      for (int i = range.last; i != last; ++ i)
	unaligned[id].insert(i);
    }
  }
  
  span_set_type spans;
  span_set_type complements;
  span_set_type unaligned;

  admissible_set_type admissibles;

  derivation_graph_type derivations;
};
