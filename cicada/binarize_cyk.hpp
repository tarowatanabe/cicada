// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_CYK__HPP__
#define __CICADA__BINARIZE_CYK__HPP__ 1

#include <deque>

#include <cicada/binarize_base.hpp>
#include <cicada/span_node.hpp>

namespace cicada
{
  struct BinarizeCYK : public BinarizeBase
  {
    BinarizeCYK(const int __order=0)
      : order(__order) {}
    
    typedef std::vector<bool, std::allocator<bool> > visited_type;
    
    typedef std::pair<int, int> span_type;
    typedef std::vector<span_type, std::allocator<span_type> > span_set_type;
    
    typedef std::vector<std::string, std::allocator<std::string> > label_set_type;
    typedef std::deque<label_set_type, std::allocator<label_set_type> > label_nodes_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > category_set_type;
    
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > ancestor_set_type;
    typedef std::deque<ancestor_set_type, std::allocator<ancestor_set_type> > ancestor_nodes_type;
    
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> >  node_set_type;
    typedef utils::chart<node_set_type, std::allocator<node_set_type> > node_chart_type;

    typedef std::vector<int, std::allocator<int> > middle_set_type;
    typedef std::deque<middle_set_type, std::allocator<middle_set_type> > middle_nodes_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, topological-order traversal to compute span for each node and 
      // ancestors.
      
      // unary rule: keep unary rules, and use the bottom as our binarized label.
      // terminal rule: keep terminal rules.

      // how to propagate features originally in the source hypergraph?
      
      // order <= 0 implies all the ancestors...

      // initialize target-nodes...
      target.clear();
      
      if (! source.is_valid()) return;
      
      // annotate spans...
      spans.clear();
      spans.reserve(source.nodes.size());
      spans.resize(source.nodes.size(), span_type(0, 0));
      
      span_node(source, spans);
      
      // start traversing in bottom-up fashion...
      labels.clear();
      labels.resize(source.nodes.size());
      
      ancestors.clear();
      ancestors.resize(source.nodes.size());

      parents.clear();
      parents.resize(source.nodes.size(), hypergraph_type::invalid);

      middles.clear();
      middles.resize(source.nodes.size());
      
      nodes.clear();
      nodes.reserve(spans.back().second + 1);
      nodes.resize(spans.back().second + 1);

      terminals.clear();
      terminals.reserve(spans.back().second);
      terminals.resize(spans.back().second);
      
      target.nodes.resize(source.nodes.size());
      target.goal = source.goal;
      for (size_t i = 0; i != target.nodes.size(); ++ i) {
	target.nodes[i].id = i;
	
	const span_type& span = spans[i];
	nodes(span.first, span.second).push_back(i);
	
	const hypergraph_type::node_type& node = source.nodes[i];
	if (! node.edges.empty()) {
	  labels[i] = label_set_type(1, source.edges[node.edges.front()].rule->lhs.non_terminal_strip());
	  
	  if (span.first + 1 == span.second)
	    terminals[span.first] = labels[i].front();
	}
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = source.edges[*eiter];
	  
	  if (edge.tails.size() <= 1) {
	    hypergraph_type::edge_type& edge_new = target.add_edge(edge.tails.begin(), edge.tails.end());
	    edge_new.rule       = edge.rule;
	    edge_new.features   = edge.features;
	    edge_new.attributes = edge.attributes;
	    
	    target.connect_edge(edge_new.id, i);
	  }
	  
	  hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	  for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	    parents[*titer] = i;
	}
      }
      
      // collect ancestors...
      if (order <= 0) {
	for (size_t i = 0; i != target.nodes.size(); ++ i) {
	  hypergraph_type::id_type parent = parents[i];
	  
	  while (parent != hypergraph_type::invalid) {
	    ancestors[i].push_back(parent);
	    parent = parents[parent];
	  }
	}
      } else {
	for (size_t i = 0; i != target.nodes.size(); ++ i) {
	  hypergraph_type::id_type parent = parents[i];

	  for (int j = 0; j != order && parent != hypergraph_type::invalid; ++ j) {
	    ancestors[i].push_back(parent);
	    parent = parents[parent];
	  }
	}
      }
      
      
      // now enumerate cyk chart...
      const int length = spans.back().second;
      
      for (int k = 2; k <= length; ++ k)
	for (int first = 0; first + k <= length; ++ first) {
	  const int last = first + k;
	  
	  for (int middle = first + 1; middle != last; ++ middle) {
	    
	    if (nodes(first, middle).empty() || nodes(middle, last).empty()) continue;
	    
	    const hypergraph_type::id_type left   = nodes(first, middle).back();
	    const hypergraph_type::id_type right  = nodes(middle, last).back();
	    
	    intersected.clear();
	    std::set_intersection(ancestors[left].begin(), ancestors[left].end(), ancestors[right].begin(), ancestors[right].end(), std::back_inserter(intersected));
	    
	    if (intersected.empty()) continue;
	    
	    if (nodes(first, last).empty()) {
	      nodes(first, last).push_back(target.add_node().id);
	      
	      ancestors.push_back(ancestor_set_type());
	      labels.push_back(label_set_type(terminals.begin() + first, terminals.begin() + last));
	      middles.push_back(middle_set_type());
	    }
	    
	    const hypergraph_type::id_type parent = nodes(first, last).front();
	    
	    unioned.clear();
	    std::set_union(ancestors[parent].begin(), ancestors[parent].end(), intersected.begin(), intersected.end(), std::back_inserter(unioned));
	    ancestors[parent].swap(unioned);
	    
	    if (labels[parent].size() > labels[left].size() + labels[right].size()) {
	      labels[parent] = labels[left];
	      labels[parent].insert(labels[parent].end(), labels[right].begin(), labels[right].end());
	    }
	    
	    middles[parent].push_back(middle);
	  }
	}
      
      // finally, register hyperedges into target!

      categories.clear();
      categories.resize(labels.size());
      
      hypergraph_type::edge_type::node_set_type tails(2);
      rule_type::symbol_set_type rhs(2);
      for (int k = 2; k <= length; ++ k)
	for (int first = 0; first + k <= length; ++ first) {
	  const int last = first + k;
	  
	  if (nodes(first, last).empty()) continue;
	  
	  const hypergraph_type::id_type parent = nodes(first, last).front();
	  
	  const symbol_type lhs = category(parent);
	  
	  middle_set_type::const_iterator miter_end = middles[parent].end();
	  for (middle_set_type::const_iterator miter = middles[parent].begin(); miter != miter_end; ++ miter) {
	    const int& middle = *miter;
	    
	    if (nodes(first, middle).empty() || nodes(middle, last).empty()) continue;
	    
	    const hypergraph_type::id_type left   = nodes(first, middle).back();
	    const hypergraph_type::id_type right  = nodes(middle, last).back();
	    
	    tails.front() = left;
	    tails.back()  = right;
	    
	    rhs.front() = category(left);
	    rhs.back()  = category(right);
	    
	    hypergraph_type::edge_type& edge = target.add_edge(tails.begin(), tails.end());
	    edge.rule = rule_type::create(rule_type(lhs, rhs.begin(), rhs.end()));
	    
	    // what features assigned???
	    
	    target.connect_edge(edge.id, parent);
	  }
	}

      hypergraph_type sorted;
      if (target.is_valid())
	topologically_sort(target, sorted);
      target.swap(sorted);
    }

    const symbol_type& category(hypergraph_type::id_type id)
    {
      if (categories[id].empty()) {
	if (! labels[id].empty()) {

	  if (labels[id].size() == 1)
	    categories[id] = '[' + labels[id].front() + ']';
	  else {
	    std::ostringstream stream;
	    stream << '[';
	    std::copy(labels[id].begin(), labels[id].end() - 1, std::ostream_iterator<std::string>(stream, "+"));
	    stream << labels[id].back() << '^'; // add indicator for binarization
	    stream << ']';
	    categories[id] = stream.str();
	  }
	} else
	  categories[id] = "[x]"; // this should not happen...
      }
      return categories[id];
    }
    
    span_set_type spans;
    
    label_set_type      terminals;
    label_nodes_type    labels;

    category_set_type   categories;
    
    ancestor_set_type   parents;
    ancestor_set_type   intersected;
    ancestor_set_type   unioned;
    ancestor_nodes_type ancestors;

    
    middle_nodes_type   middles;

    node_chart_type     nodes;
    
    const int order;
  };
  
  inline
  void binarize_cyk(const HyperGraph& source, HyperGraph& target, const int order=0)
  {
    BinarizeCYK binarizer(order);
    
    binarizer(source, target);
  }

  inline
  void binarize_cyk(HyperGraph& source, const int order=0)
  {
    HyperGraph target;

    BinarizeCYK binarizer(order);
    
    binarizer(source, target);
    
    source.swap(target);
  }
};

#endif
