// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_ALL__HPP__
#define __CICADA__BINARIZE_ALL__HPP__ 1

#include <cicada/binarize_base.hpp>

#include <utils/sgi_hash_map.hpp>
#include <utils/hashmurmur.hpp>

namespace cicada
{
  struct BinarizeAll : public BinarizeBase
  {
    
    //
    // sub-tails to label_sequence + label mapping
    // sub-tails + sub-symbols to node-id mapping
    //

    typedef rule_type::symbol_set_type                 symbol_set_type;
    typedef hypergraph_type::edge_type::node_set_type  tail_set_type;
    typedef std::pair<tail_set_type, symbol_set_type > tail_symbol_pair_type;

    struct tail_set_hash : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      size_t operator()(const tail_set_type& x) const {
	return hasher_type::operator()(x.begin(), x.end(), 0);
      }
    };

    struct tail_symbol_pair_hash : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const tail_symbol_pair_type& x) const {
	return hasher_type::operator()(x.first.begin(), x.first.end(), hasher_type::operator()(x.second.begin(), x.second.end(), 0));
      }
    };

#ifdef HAVE_TR1_UNORDERED_MAP
      typedef std::tr1::unordered_map<tail_set_type, symbol_type, tail_set_hash, std::equal_to<tail_set_type>,
				      std::allocator<std::pair<const tail_set_type, symbol_type> > > label_map_type;
      typedef std::tr1::unordered_map<tail_symbol_pair_type, hypergraph_type::id_type, tail_symbol_pair_hash, std::equal_to<tail_symbol_pair_type>,
				      std::allocator<std::pair<const tail_symbol_pair_type, hypergraph_type::id_type> > > node_map_type;
#else
      typedef sgi::hash_map<tail_set_type, symbol_type, tail_set_hash, std::equal_to<tail_set_type>,
			    std::allocator<std::pair<const tail_set_type, symbol_type> > > label_map_type;
      typedef sgi::hash_map<tail_symbol_pair_type, hypergraph_type::id_type, tail_symbol_pair_hash, std::equal_to<tail_symbol_pair_type>,
			    std::allocator<std::pair<const tail_symbol_pair_type, hypergraph_type::id_type> > > node_map_type;
#endif

    label_map_type label_map;
    node_map_type node_map;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, copy...
      target = source;

      if (! source.is_valid()) return;

      phrase_type       binarized(2);
      hypergraph_type::edge_type::node_set_type tails(2);
      
      // we will traverse source-side in order to avoid confusion with newly created nodes...
      removed_type removed(source.edges.size(), false);

      position_set_type positions;
      node_chart_type   node_chart;
      label_chart_type  label_chart;
      
      label_map.clear();
      node_map.clear();
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  
	  if (edge_source.tails.size() <= 2) continue;
	  
	  removed[edge_source.id] = true;
	  
	  // we will create nodes in a chart structure, and exhaustively enumerate edges
	  
	  symbol_set_type rhs_sorted(edge_source.rule->rhs);
	  tail_set_type   tails_sorted(edge_source.tails);
	  
	  // first, compute non-terminal spans...
	  positions.clear();
	  int pos = 0;
	  for (size_t i = 0; i != rhs_sorted.size(); ++ i)
	    if (rhs_sorted[i].is_non_terminal()) {
	      const int non_terminal_index = rhs_sorted[i].non_terminal_index();
	      
	      tails_sorted[pos] = edge_source.tails[utils::bithack::branch(non_terminal_index == 0, pos, non_terminal_index - 1)];
	      
	      rhs_sorted[i] = rhs_sorted[i].non_terminal();
	      
	      positions.push_back(i);
	      
	      ++ pos;
	    }
	  
	  if (positions.size() != edge_source.tails.size())
	    throw std::runtime_error("invalid edge: # of non-terminals and tails size do not match");
	  
	  // seond, enumerate chart to compute node and edges...
	  node_chart.clear();
	  node_chart.resize(positions.size() + 1, hypergraph_type::invalid);

	  label_chart.clear();
	  label_chart.resize(positions.size() + 1);
	  
	  for (size_t i = 0; i != positions.size(); ++ i) {
	    node_chart(i, i + 1) = tails_sorted[i];
	    label_chart(i, i + 1) = rhs_sorted[positions[i]];
	  }
	  
	  for (size_t length = 2; length < positions.size(); ++ length)
	    for (size_t first = 0; first + length <= positions.size(); ++ first) {
	      const size_t last = first + length;
	      
	      const symbol_set_type subrhs(rhs_sorted.begin() + positions[first], rhs_sorted.begin() + positions[last - 1] + 1);
	      const tail_set_type   subtails(tails_sorted.begin() + first, tails_sorted.begin() + last);
	      
	      std::pair<label_map_type::iterator, bool> result_label = label_map.insert(std::make_pair(subtails, symbol_type()));
	      if (result_label.second) {
		const symbol_type::piece_type left = label_chart(first, last - 1).non_terminal_strip();
		const symbol_type::piece_type right = label_chart(last - 1, last).non_terminal_strip();

		if (length > 2)
		  result_label.first->second = '[' + std::string(left.begin(), left.end() - 1) + '+' + right + "^]";
		else
		  result_label.first->second = '[' + std::string(left) + '+' + right + "^]";
	      }
	      
	      std::pair<node_map_type::iterator, bool> result_node = node_map.insert(std::make_pair(tail_symbol_pair_type(subtails, subrhs), 0));
	      if (result_node.second)
		result_node.first->second = target.add_node().id;
	      
	      const symbol_type lhs = result_label.first->second;
	      const hypergraph_type::id_type head = result_node.first->second;
	      
	      node_chart(first, last) = head;
	      label_chart(first, last) = lhs;
	      
	      // if newly created, then, create edges
	      if (result_node.second)
		for (size_t middle = first + 1; middle != last; ++ middle) {
		  // [first, middle) and [middle, last)
		  
		  tails.front() = node_chart(first, middle);
		  tails.back()  = node_chart(middle, last);
		  
		  const size_t middle_first = positions[middle - 1] + 1;
		  const size_t middle_last  = positions[middle];
		  
		  binarized.clear();
		  binarized.push_back(label_chart(first, middle));
		  binarized.insert(binarized.end(), rhs_sorted.begin() + middle_first, rhs_sorted.begin() + middle_last);
		  binarized.push_back(label_chart(middle, last));
		  
		  hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
		  edge_new.rule = rule_type::create(rule_type(lhs, binarized.begin(), binarized.end()));
		  
		  target.connect_edge(edge_new.id, head);
		}
	    }
	  
	  // root...
	  {
	    const size_t first = 0;
	    const size_t last = positions.size();
	    
	    const hypergraph_type::id_type head = node_source.id;
	    const symbol_type& lhs = edge_source.rule->lhs;
	    
	    node_chart(first, last) = head;
	    label_chart(first, last) = lhs;
	    
	    for (size_t middle = first + 1; middle != last; ++ middle) {
	      // [first, middle) and [middle, last)
	      
	      tails.front() = node_chart(first, middle);
	      tails.back()  = node_chart(middle, last);
	      
	      binarized.clear();
	      
	      const size_t prefix_first = 0;
	      const size_t prefix_last  = positions[first];
	      
	      binarized.insert(binarized.end(), rhs_sorted.begin() + prefix_first, rhs_sorted.begin() + prefix_last);
	      binarized.push_back(label_chart(first, middle));
	      
	      const size_t middle_first = positions[middle - 1] + 1;
	      const size_t middle_last  = positions[middle];
	      
	      binarized.insert(binarized.end(), rhs_sorted.begin() + middle_first, rhs_sorted.begin() + middle_last);
	      binarized.push_back(label_chart(middle, last));
	      
	      const size_t suffix_first = positions[last - 1] + 1;
	      const size_t suffix_last  = rhs_sorted.size();
	      
	      binarized.insert(binarized.end(), rhs_sorted.begin() + suffix_first, rhs_sorted.begin() + suffix_last);
	      
	      hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	      edge_new.rule       = rule_type::create(rule_type(lhs, binarized.begin(), binarized.end()));
	      edge_new.features   = edge_source.features;
	      edge_new.attributes = edge_source.attributes;
	      
	      target.connect_edge(edge_new.id, head);
	    }
	  }
	}
      }
      
      // further resize...
      removed.resize(target.edges.size(), false);
      
      hypergraph_type graph_removed;
      
      topologically_sort(target, graph_removed, filter(removed));
      
      target.swap(graph_removed);
    }
  };

  inline
  void binarize_all(const HyperGraph& source, HyperGraph& target)
  {
    BinarizeAll binarizer;
    
    binarizer(source, target);
  }

  inline
  void binarize_all(HyperGraph& source)
  {
    HyperGraph target;

    BinarizeAll binarizer;
    
    binarizer(source, target);
    
    source.swap(target);
  }  
};

#endif
