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

    template <typename Iterator>
    symbol_type binarized_label(const symbol_type& lhs, Iterator first, Iterator last)
    {
      std::string label = lhs.non_terminal().non_terminal_strip() + '^';
      for (/**/; first != last; ++ first) {
	label += '_';
	label += first->non_terminal().non_terminal_strip();
      }
      return '[' + label + ']';
    }
    
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
	  
	  removed[edge_source.id] = true;
	  
	  // right most antecedents binarization...
	  
	  binarized.clear();
	  context.clear();
	  
	  hypergraph_type::id_type head = edge_source.head;
	  symbol_type non_terminal_head = edge_source.rule->lhs;
	  
	  const int arity = edge_source.tails.size();
	  
	  int pos = 0;
	  rule_type::symbol_set_type::const_iterator riter = edge_source.rule->rhs.begin();
	  rule_type::symbol_set_type::const_iterator riter_end = edge_source.rule->rhs.end();
	  for (/**/; riter != riter_end && pos < arity - 2; ++ riter) {
	    if (riter->is_non_terminal()) {
	      // we will check non-terminal-index here...
	      const int non_terminal_index = riter->non_terminal_index();
	      if (non_terminal_index > 0 && non_terminal_index - 1 != pos)
		throw std::runtime_error("hypergraph is not tail-sorted!");
	      
	      // then...
	      context.push_back(*riter);
	      
	      const symbol_type non_terminal_new = (order < 0
						    ? binarized_label(edge_source.rule->lhs, context.begin(), context.end())
						    : binarized_label(edge_source.rule->lhs,
								      std::max(context.begin(), context.end() - order),
								      context.end()));
	      
	      binarized.push_back(riter->non_terminal());
	      binarized.push_back(non_terminal_new);
	      
	      hypergraph_type::node_type& node_new = target.add_node();
	      tails.front() = edge_source.tails[pos];
	      tails.back() = node_new.id;
	      
	      hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	      
	      edge_new.rule = rule_type::create(rule_type(non_terminal_head, binarized.begin(), binarized.end()));
	      
	      target.connect_edge(edge_new.id, head);
	      
	      head = node_new.id;
	      non_terminal_head = non_terminal_new;
	      
	      binarized.clear();
	      ++ pos;
	    } else 
	      binarized.push_back(*riter);
	  }
	  
	  hypergraph_type::edge_type& edge_new = target.add_edge(edge_source.tails.end() - 2, edge_source.tails.end());
	  
	  for (/**/; riter != riter_end; ++ riter) {
	    if (riter->is_non_terminal()) {
	      // we will check non-terminal-index here...
	      const int non_terminal_index = riter->non_terminal_index();
	      if (non_terminal_index > 0 && non_terminal_index - 1 != pos)
		throw std::runtime_error("hypergraph is not tail-sorted!");
	      
	      ++ pos;
	      binarized.push_back(riter->non_terminal());
	    } else
	      binarized.push_back(*riter);
	  }
	  
	  edge_new.rule = rule_type::create(rule_type(non_terminal_head, binarized.begin(), binarized.end()));
	  
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
    
    template <typename Iterator>
    symbol_type binarized_label(const symbol_type& lhs, Iterator first, Iterator last)
    {
      std::string label = lhs.non_terminal().non_terminal_strip() + '^';
      for (/**/; first != last; ++ first) {
	label += '~';
	label += first->non_terminal().non_terminal_strip();
      }
      return '[' + label + ']';
    }
    
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
	  
	  removed[edge_source.id] = true;
	  
	  // left most antecedents binarization...
	  
	  binarized.clear();
	  context.clear();
	  
	  hypergraph_type::id_type head = edge_source.head;
	  symbol_type non_terminal_head = edge_source.rule->lhs;
	  
	  const int arity = edge_source.tails.size();

	  int pos = arity - 1;
	  rule_type::symbol_set_type::const_reverse_iterator riter = edge_source.rule->rhs.rbegin();
	  rule_type::symbol_set_type::const_reverse_iterator riter_end = edge_source.rule->rhs.rend();
	  for (/**/; riter != riter_end && pos >= 2; ++ riter) {
	    if (riter->is_non_terminal()) {
	      // we will check non-terminal-index here...
	      const int non_terminal_index = riter->non_terminal_index();
	      if (non_terminal_index > 0 && non_terminal_index - 1 != pos)
		throw std::runtime_error("hypergraph is not tail-sorted!");
	      
	      context.push_back(*riter);
	      
	      const symbol_type non_terminal_new = (order < 0
						    ? binarized_label(edge_source.rule->lhs, context.begin(), context.end())
						    : binarized_label(edge_source.rule->lhs,
								      std::max(context.begin(), context.end() - order),
								      context.end()));
	      
	      binarized.push_back(*riter);
	      binarized.push_back(non_terminal_new);
	      
	      hypergraph_type::node_type& node_new = target.add_node();
	      tails.front() =  node_new.id;
	      tails.back() = edge_source.tails[pos];
	      
	      hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());

	      std::reverse(binarized.begin(), binarized.end());
	      edge_new.rule = rule_type::create(rule_type(non_terminal_head, binarized.begin(), binarized.end()));
	      
	      target.connect_edge(edge_new.id, head);
	      
	      head = node_new.id;
	      non_terminal_head = non_terminal_new;
	      
	      binarized.clear();
	      -- pos;
	    } else
	      binarized.push_back(*riter);
	  }
	  
	  hypergraph_type::edge_type& edge_new = target.add_edge(edge_source.tails.begin(), edge_source.tails.begin() + 2);

	  for (/**/; riter != riter_end; ++ riter) {
	    if (riter->is_non_terminal()) {
	      // we will check non-terminal-index here...
	      const int non_terminal_index = riter->non_terminal_index();
	      if (non_terminal_index > 0 && non_terminal_index - 1 != pos)
		throw std::runtime_error("hypergraph is not tail-sorted!");
	      
	      -- pos;
	      binarized.push_back(riter->non_terminal());
	    } else
	      binarized.push_back(*riter);
	  }
	  
	  std::reverse(binarized.begin(), binarized.end());
	  
	  edge_new.rule = rule_type::create(rule_type(non_terminal_head, binarized.begin(), binarized.end()));
	  
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
