// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_TERMINAL__HPP__
#define __CICADA__BINARIZE_TERMINAL__HPP__ 1

#include <cicada/binarize_base.hpp>

namespace cicada
{
  
  struct BinarizeTerminal : public BinarizeBase
  {
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, copy...
      target = source;
      
      phrase_type binarized;
      tails_type  tails;
      
      phrase_type binarized_local;
      tails_type  tails_local;
      
      // we will traverse source-side in order to avoid confusion with newly created nodes...
      removed_type removed(source.edges.size(), false);
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  
	  if (edge_source.tails.size() <= 1) continue;
	  	  
	  // we will create nodes in a chart structure, and exhaustively enumerate edges
	  
	  const rule_type::symbol_set_type& rhs = edge_source.rule->rhs;
	  const symbol_type& lhs = edge_source.rule->lhs;
	  const symbol_type lhs_binarized = '[' + edge_source.rule->lhs.non_terminal_strip() + "^]";
	  
	  binarized.clear();
	  tails.clear();
	  tails_local.clear();
	  int pos = 0;
	  for (size_t i = 0; i != rhs.size(); ++ i)
	    if (rhs[i].is_non_terminal()) {
	      const symbol_type& non_terminal = rhs[i];
	      
	      // we will check non-terminal-index here...
	      const int non_terminal_index = non_terminal.non_terminal_index();
	      if (non_terminal_index > 0 && non_terminal_index - 1 != pos)
		throw std::runtime_error("hypergraph is not tail-sorted!");
	      
	      tails_local.push_back(edge_source.tails[pos]);
	      
	      ++ pos;
	    } else {
	      if (! tails_local.empty()) {
		if (tails_local.size() == 1) {
		  tails.push_back(tails_local.front());
		  binarized.push_back(rhs[i - 1].non_terminal());
		} else {
		  const hypergraph_type::id_type node_id = binarize(lhs_binarized, lhs_binarized, hypergraph_type::invalid,
								    feature_set_type(), attribute_set_type(),
								    tails_local.begin(), tails_local.end(),
								    rhs.begin() + i - tails_local.size(), rhs.begin() + i,
								    target);
		  
		  tails.push_back(node_id);
		  binarized.push_back(lhs_binarized);
		}
		tails_local.clear();
	      }
	      binarized.push_back(rhs[i]);
	    }
	  
	  if (binarized.empty()) {
	    // this will happen when no terminals found...
	    
	    // we will not binarize when tails_local.size() == 2...
	    if (tails_local.size() > 2) {
	      binarize(lhs, lhs_binarized, node_source.id,
		       edge_source.features, edge_source.attributes,
		       tails_local.begin(), tails_local.end(),
		       rhs.end() - tails_local.size(), rhs.end(),
		       target);
	      
	      removed[edge_source.id] = true;
	    }
	  } else {
	    if (! tails_local.empty()) {
	      if (tails_local.size() == 1) {
		tails.push_back(tails_local.back());
		binarized.push_back(rhs.back());
	      } else {
		const hypergraph_type::id_type node_id = binarize(lhs_binarized, lhs_binarized, hypergraph_type::invalid,
								  feature_set_type(), attribute_set_type(),
								  tails_local.begin(), tails_local.end(),
								  rhs.end() - tails_local.size(), rhs.end(),
								  target);
		
		tails.push_back(node_id);
		binarized.push_back(lhs_binarized);
	      }
	    }
	    
	    // create an edge leadning to node_soucrce.id...
	    hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	    
	    edge_new.rule = rule_type::create(rule_type(lhs, binarized.begin(), binarized.end()));
	    edge_new.features   = edge_source.features;
	    edge_new.attributes = edge_source.attributes;
	    
	    target.connect_edge(edge_new.id, node_source.id);
	    
	    removed[edge_source.id] = true;
	  }
	}
      }
      
      // further resize...
      removed.resize(target.edges.size(), false);
      
      hypergraph_type graph_removed;
      
      topologically_sort(target, graph_removed, filter(removed));
      
      target.swap(graph_removed);
    }
    
  private:
    
    node_chart_type   node_chart;
    label_chart_type  label_chart;
    
    template <typename IteratorTail, typename IteratorPhrase>
    hypergraph_type::id_type binarize(const symbol_type& lhs,
				      const symbol_type& lhs_binarized,
				      const hypergraph_type::id_type& head,
				      const feature_set_type& features,
				      const attribute_set_type& attributes,
				      IteratorTail tail_first, IteratorTail tail_last,
				      IteratorPhrase phrase_first, IteratorPhrase phrase_last,
				      hypergraph_type& graph)
    {
      const size_t tail_size = tail_last - tail_first;

      node_chart.clear();
      node_chart.resize(tail_size + 1, hypergraph_type::invalid);
      
      label_chart.clear();
      label_chart.resize(tail_size + 1);
      
      for (size_t i = 0; i != tail_size; ++ i) {
	node_chart(i, i + 1) = tail_first[i];
	label_chart(i, i + 1) = phrase_first[i].non_terminal();
      }
      
      phrase_type       binarized(2);
      hypergraph_type::edge_type::node_set_type tails(2);
      
      for (size_t length = 2; length <= tail_size; ++ length)
	for (size_t first = 0; first + length <= tail_size; ++ first) {
	  const bool is_root = length == tail_size;
	  const size_t last = first + length;
	  
	  const hypergraph_type::id_type head_rule = (is_root && head != hypergraph_type::invalid ? head : graph.add_node().id);
	  const symbol_type& lhs_rule = (is_root ? lhs : lhs_binarized);
	  
	  node_chart(first, last) = head_rule;
	  label_chart(first, last) = lhs_rule;
	  
	  for (size_t middle = first + 1; middle != last; ++ middle) {
	    // [first, middle) and [middle, last)
	    
	    tails.front() = node_chart(first, middle);
	    tails.back()  = node_chart(middle, last);
	    
	    binarized.front() = label_chart(first, middle);
	    binarized.back() = label_chart(middle, last);
	    
	    hypergraph_type::edge_type& edge_new = graph.add_edge(tails.begin(), tails.end());
	    edge_new.rule = rule_type::create(rule_type(lhs, binarized.begin(), binarized.end()));
	    
	    graph.connect_edge(edge_new.id, head_rule);
	    
	    // assign features...
	    if (is_root) {
	      edge_new.features   = features;
	      edge_new.attributes = attributes;
	    }
	  }
	}
      
      return node_chart(0, tail_size);
    }
  };

  inline
  void binarize_terminal(const HyperGraph& source, HyperGraph& target)
  {
    BinarizeTerminal binarizer;
    
    binarizer(source, target);
  }

  inline
  void binarize_terminal(HyperGraph& source)
  {
    HyperGraph target;

    BinarizeTerminal binarizer;
    
    binarizer(source, target);
    
    source.swap(target);
  }

};

#endif
