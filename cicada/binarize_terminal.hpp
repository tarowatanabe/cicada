// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_TERMINAL__HPP__
#define __CICADA__BINARIZE_TERMINAL__HPP__ 1

#include <cicada/binarize_base.hpp>

// we assume that every hyperedge takes the structure of:
// 
// [x]([x]^* terminal [x]^*)
//
// it is binarized so that:
// [x]([x^L] terminal [x^R])
// becomes
// [x^L]([x^L][x]) // left heavy...
// [x^L]([x] [x])
// [x^R]([x][x^R]) // right heavy...
// [x^R]([x] [x])

namespace cicada
{
  
  struct BinarizeTerminal : public BinarizeBase
  {
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > rhs_type;

    typedef rule_type::symbol_set_type                 symbol_set_type;
    typedef hypergraph_type::edge_type::node_set_type  tail_set_type;

    typedef std::pair<symbol_type, tail_set_type> binarized_type;
    
    struct binarized_hash : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const binarized_type& x) const {
	return hasher_type::operator()(x.second.begin(), x.second.end(), x.first.id());
      }
    };
    
    typedef utils::unordered_map<binarized_type, hypergraph_type::id_type, binarized_hash, std::equal_to<binarized_type>,
				 std::allocator<std::pair<const binarized_type, hypergraph_type::id_type> > >::type binarized_map_type;

    binarized_map_type binarized_left;
    binarized_map_type binarized_right;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, copy...
      target = source;
      
      if (! source.is_valid()) return;
      
      symbol_type terminal;
      rhs_type    rhs_left;
      rhs_type    rhs_right;
      tails_type  tails_left;
      tails_type  tails_right;

      rhs_type   rhs_new;
      tails_type tails_new;
      
      removed_type removed(source.edges.size(), false);

      binarized_left.clear();
      binarized_right.clear();
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  
	  if (edge_source.tails.size() <= 1) continue;
	  
	  terminal = symbol_type();
	  rhs_left.clear();
	  rhs_right.clear();
	  tails_left.clear();
	  tails_right.clear();
	  
	  int non_terminal_pos = 0;
	  
	  rule_type::symbol_set_type::const_iterator riter_end = edge_source.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator riter = edge_source.rule->rhs.begin(); riter != riter_end; ++ riter)
	    if (riter->is_non_terminal()) {
	      const int __non_terminal_index = riter->non_terminal_index();
	      const int non_terminal_index = utils::bithack::branch(__non_terminal_index <= 0,
								    non_terminal_pos,
								    __non_terminal_index - 1);
	      ++ non_terminal_pos;
	      
	      if (terminal == symbol_type()) {
		rhs_left.push_back(riter->non_terminal());
		tails_left.push_back(edge_source.tails[non_terminal_index]);
	      } else {
		rhs_right.push_back(riter->non_terminal());
		tails_right.push_back(edge_source.tails[non_terminal_index]);
	      }
	    } else {
	      if (terminal != symbol_type())
		throw std::runtime_error("invalid dependency structure");
	      
	      terminal = *riter;
	    }
	  
	  // check if we need to binarize...
	  if (tails_left.size() <= 1 && tails_right.size() <= 1) continue;

	  const symbol_type& lhs = edge_source.rule->lhs;
	  
	  rhs_new.clear();
	  tails_new.clear();
	  
	  if (! tails_left.empty()) {
	    if (tails_left.size() == 1) {
	      rhs_new.push_back(rhs_left.front());
	      tails_new.push_back(tails_left.front());
	    } else {
	      const std::pair<symbol_type, hypergraph_type::id_type> result = binarize_left(lhs, rhs_left, tails_left, target);
	      
	      rhs_new.push_back(result.first);
	      tails_new.push_back(result.second);
	    }
	  }
	  
	  rhs_new.push_back(terminal);
	  
	  if (! tails_right.empty()) {
	    if (tails_right.size() == 1) {
	      rhs_new.push_back(rhs_right.front());
	      tails_new.push_back(tails_right.front());
	    } else {
	      const std::pair<symbol_type, hypergraph_type::id_type> result = binarize_right(lhs, rhs_right, tails_right, target);
	      
	      rhs_new.push_back(result.first);
	      tails_new.push_back(result.second);
	    }
	  }
	  
	  hypergraph_type::edge_type& edge_new = target.add_edge(tails_new.begin(), tails_new.end());
	  edge_new.rule       = rule_type::create(rule_type(lhs, rhs_new.begin(), rhs_new.end()));
	  edge_new.features   = edge_source.features;
	  edge_new.attributes = edge_source.attributes;
	  
	  target.connect_edge(edge_new.id, node_source.id);
	  
	  removed[edge_source.id] = true;
	}
      }
      
      // further resize...
      removed.resize(target.edges.size(), false);
      
      hypergraph_type graph_removed;
      
      topologically_sort(target, graph_removed, filter(removed));
      
      target.swap(graph_removed);
    }
    
    std::pair<symbol_type, hypergraph_type::id_type> binarize_left(const symbol_type& lhs,
								   const rhs_type& rhs,
								   const tails_type& tails,
								   hypergraph_type& target)
    {
      // proceed in a left-to-right order, and construct...
      const symbol_type lhs_left = "[" + lhs.non_terminal_strip() + "^L]";
      
      symbol_set_type rhs_binarized(2, rhs.front());
      tail_set_type   tails_binarized(2, tails.front());
      
      rhs_type::const_iterator riter_end = rhs.end();
      tails_type::const_iterator titer = tails.begin() + 1;
      for (rhs_type::const_iterator riter = rhs.begin() + 1; riter != riter_end; ++ riter, ++ titer) {
	rhs_binarized.back() = *riter;
	tails_binarized.back() = *titer;
	
	std::pair<binarized_map_type::iterator, bool> result = binarized_left.insert(std::make_pair(binarized_type(lhs_left,
														   tails_binarized),
												  0));
	if (result.second) {
	  const hypergraph_type::id_type node_id = target.add_node().id;
	  
	  hypergraph_type::edge_type& edge_new = target.add_edge(tails_binarized.begin(), tails_binarized.end());
	  
	  edge_new.rule = rule_type::create(rule_type(lhs_left, rhs_binarized.begin(), rhs_binarized.end()));
	  
	  target.connect_edge(edge_new.id, node_id);
	  
	  result.first->second = edge_new.id;
	}
	
	rhs_binarized.front() = lhs_left;
	tails_binarized.front() = result.first->second;
      }
      
      return std::make_pair(rhs_binarized.front(), tails_binarized.front());
    }
    
    std::pair<symbol_type, hypergraph_type::id_type> binarize_right(const symbol_type& lhs,
								    const rhs_type& rhs,
								    const tails_type& tails,
								    hypergraph_type& target)
    {
      const symbol_type lhs_right = "[" + lhs.non_terminal_strip() + "^R]";
      
      symbol_set_type rhs_binarized(2, rhs.back());
      tail_set_type   tails_binarized(2, tails.back());
      
      rhs_type::const_reverse_iterator riter_end = rhs.rend();
      tails_type::const_reverse_iterator titer = tails.rbegin() + 1;
      for (rhs_type::const_reverse_iterator riter = rhs.rbegin() + 1; riter != riter_end; ++ riter, ++ titer) {
	rhs_binarized.front() = *riter;
	tails_binarized.front() = *titer;
	
	std::pair<binarized_map_type::iterator, bool> result = binarized_right.insert(std::make_pair(binarized_type(lhs_right,
														    tails_binarized),
												   0));
	
	if (result.second) {
	  const hypergraph_type::id_type node_id = target.add_node().id;
	  
	  hypergraph_type::edge_type& edge_new = target.add_edge(tails_binarized.begin(), tails_binarized.end());
	  
	  edge_new.rule = rule_type::create(rule_type(lhs_right, rhs_binarized.begin(), rhs_binarized.end()));
	  
	  target.connect_edge(edge_new.id, node_id);
	  
	  result.first->second = edge_new.id;
	}
	
	rhs_binarized.back() = lhs_right;
	tails_binarized.back() = result.first->second;
      }
      
      return std::make_pair(rhs_binarized.back(), tails_binarized.back());
    }
    
#if 0
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, copy...
      target = source;
      
      if (! source.is_valid()) return;
      
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
#endif
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
