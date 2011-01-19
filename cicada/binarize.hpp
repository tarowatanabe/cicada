// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE__HPP__
#define __CICADA__BINARIZE__HPP__ 1

#include <stdexcept>
#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/sort.hpp>

#include <utils/bithack.hpp>
#include <utils/hashmurmur.hpp>

#include <boost/lexical_cast.hpp>

#include <google/dense_hash_set>

namespace cicada
{
  
  struct __BinarizeBase
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::symbol_type      symbol_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > phrase_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > context_type;
    
    typedef std::vector<bool, std::allocator<bool> > removed_type;

    template <typename Iterator>
    std::string binarized_label(const symbol_type& lhs, Iterator first, Iterator last)
    {
      std::string label = lhs.non_terminal().non_terminal_strip() + '^';
      if (first != last) {
	label += first->non_terminal().non_terminal_strip();
	++ first;
	for (/**/; first != last; ++ first) {
	  label += '_';
	  label += first->non_terminal().non_terminal_strip();
	}
      }
      return label;
    }

    
    struct filter
    {
      const removed_type& removed;
      
      filter(const removed_type& __removed) : removed(__removed) {}
      
      template <typename Edge>
      bool operator()(const Edge& edge) const
      {
	return removed[edge.id];
      }
    };
  };

  struct BinarizeRight : public __BinarizeBase
  {
    BinarizeRight(const int __order=-1)
      : order(__order) {}


    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, copy...
      target = source;
      
      phrase_type binarized(2);
      hypergraph_type::edge_type::node_set_type tails(2);
      
      // we will traverse source-side in order to avoid confusion with newly created nodes...
      
      removed_type removed(source.edges.size(), false);

      context_type context;
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  	  
	  if (edge_source.tails.size() <= 2) continue;
	  
	  if (edge_source.tails.size() != static_cast<size_t>(edge_source.rule->rhs.size()))
	    throw std::runtime_error("we do not support terminal-mixed rules (aka Hiero rules)");
	  
	  removed[edge_source.id] = true;
	  
	  // right most antecedents binarization...
	  
	  context.clear();
	  
	  hypergraph_type::id_type head = edge_source.head;
	  std::string non_terminal_head = edge_source.rule->lhs.non_terminal_strip();
	  
	  const int arity = edge_source.tails.size();
	  
	  for (int i = 0; i < arity - 2; ++ i) {
	    context.push_back(edge_source.rule->rhs[i]);
	    
	    const std::string non_terminal_new = (order < 0
						  ? binarized_label(edge_source.rule->lhs, context.begin(), context.end())
						  : binarized_label(edge_source.rule->lhs,
								    std::max(context.begin(), context.end() - order),
								    context.end()));
	      
	    hypergraph_type::node_type& node_new = target.add_node();
	    tails.front() = edge_source.tails[i];
	    tails.back() = node_new.id;
	    
	    binarized.front() = edge_source.rule->rhs[i];
	    binarized.back() = '[' + non_terminal_new + ']';
	    
	    hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	    
	    edge_new.rule = rule_type::create(rule_type('[' + non_terminal_head + ']', binarized.begin(), binarized.end()));
	    
	    target.connect_edge(edge_new.id, head);
	    
	    head = node_new.id;
	    non_terminal_head = non_terminal_new;
	  }
	  
	  hypergraph_type::edge_type& edge_new = target.add_edge(edge_source.tails.end() - 2, edge_source.tails.end());
	  
	  binarized.front() = *(edge_source.rule->rhs.end() - 2);
	  binarized.back()  = *(edge_source.rule->rhs.end() - 1);
	  
	  edge_new.rule = rule_type::create(rule_type('[' + non_terminal_head + ']', binarized.begin(), binarized.end()));
	  
	  // assign features here...
	  edge_new.features   = edge_source.features;
	  edge_new.attributes = edge_source.attributes;
	  
	  target.connect_edge(edge_new.id, head);
	}
      }
      
      // further resize...
      removed.resize(target.edges.size(), false);
      
      hypergraph_type graph_removed;
      
      topologically_sort(target, graph_removed, filter(removed));
      
      target.swap(graph_removed);
    }
    
    int order;
  };


  struct BinarizeLeft : public __BinarizeBase
  {
    BinarizeLeft(const int __order=-1)
      : order(__order) {}
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, copy...
      target = source;
      
      phrase_type binarized(2);
      hypergraph_type::edge_type::node_set_type tails(2);

      context_type context;
      
      // we will traverse source-side in order to avoid confusion with newly created nodes...
      
      removed_type removed(source.edges.size(), false);
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  
	  if (edge_source.tails.size() <= 2) continue;
	  
	  if (edge_source.tails.size() != static_cast<size_t>(edge_source.rule->rhs.size()))
	    throw std::runtime_error("we do not support terminal-mixed rules (aka Hiero rules)");
	  
	  removed[edge_source.id] = true;
	  
	  // left most antecedents binarization...

	  context.clear();
	  
	  hypergraph_type::id_type head = edge_source.head;
	  std::string non_terminal_head = edge_source.rule->lhs.non_terminal_strip();
	  
	  const int arity = edge_source.tails.size();
	  
	  for (int i = 0; i < arity - 2; ++ i) {
	    context.push_back(edge_source.rule->rhs[arity - i - 1]);
	    
	    const std::string non_terminal_new = (order < 0
						  ? binarized_label(edge_source.rule->lhs, context.begin(), context.end())
						  : binarized_label(edge_source.rule->lhs,
								    std::max(context.begin(), context.end() - order),
								    context.end()));
	    
	    hypergraph_type::node_type& node_new = target.add_node();
	    tails.front() =  node_new.id;
	    tails.back() = edge_source.tails[arity - i - 1];
	    
	    binarized.front() =  '[' + non_terminal_new + ']';
	    binarized.back() = edge_source.rule->rhs[arity - i - 1];
	    
	    hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	    
	    edge_new.rule = rule_type::create(rule_type('[' + non_terminal_head + ']', binarized.begin(), binarized.end()));
	    
	    target.connect_edge(edge_new.id, head);
	    
	    head = node_new.id;
	    non_terminal_head = non_terminal_new;
	  }
	  
	  hypergraph_type::edge_type& edge_new = target.add_edge(edge_source.tails.begin(), edge_source.tails.begin() + 1);
	  
	  binarized.front() = *(edge_source.rule->rhs.begin() + 0);
	  binarized.back()  = *(edge_source.rule->rhs.begin() + 1);
	  
	  edge_new.rule = rule_type::create(rule_type('[' + non_terminal_head + ']', binarized.begin(), binarized.end()));
	  
	  // assign features here...
	  edge_new.features   = edge_source.features;
	  edge_new.attributes = edge_source.attributes;
	  
	  target.connect_edge(edge_new.id, head);
	}
      }
      
      // further resize...
      removed.resize(target.edges.size(), false);
      
      hypergraph_type graph_removed;
      
      topologically_sort(target, graph_removed, filter(removed));
      
      target.swap(graph_removed);
    }
    
    int order;
  };
  
  
  inline
  void binarize_right(const HyperGraph& source, HyperGraph& target, const int order=-1)
  {
    BinarizeRight binarizer(order);
    
    binarizer(source, target);
  }

  inline
  void binarize_right(HyperGraph& source, const int order=-1)
  {
    HyperGraph target;

    BinarizeRight binarizer(order);
    
    binarizer(source, target);
    
    source.swap(target);
  }

  inline
  void binarize_left(const HyperGraph& source, HyperGraph& target, const int order=-1)
  {
    BinarizeLeft binarizer(order);
    
    binarizer(source, target);
  }

  inline
  void binarize_left(HyperGraph& source, const int order=-1)
  {
    HyperGraph target;

    BinarizeLeft binarizer(order);
    
    binarizer(source, target);
    
    source.swap(target);
  }

};

#endif
