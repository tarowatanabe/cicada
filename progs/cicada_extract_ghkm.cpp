
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
      if (span.empty())
	return false;
      else
	return span.lower_bound(range.first) != span.lower_bound(range.last);
    }
  
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
    span_set_type spans(graph.nodes.size());
    span_set_type complements(graph.nodes.size());
    span_set_type unaligned(graph.nodes.size());
    
    admissible_set_type admissibles(graph.nodes.size());
    
    admissible_nodes(graph, sentence, alignment, spans, admissibles);

    
    minimal_rules(graph, spans, admissibles)
    
  }
  
  void minimal_rules(const hypergraph_type& graph,
		     const span_set_type& spans,
		     const admissible_set_type& admissibles)
  {
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
			const alignment_type& alignment,
			span_set_type& spans,
			span_set_type& complements,
			admissible_set_type& admissibles)
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
    for (size_t id = 0; id != graph.nodes.size(); ++ id)
      admissibles[id] = ! complements[id].intersect(spans[id].range());
  }
};
