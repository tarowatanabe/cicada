// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_LEFT__HPP__
#define __CICADA__BINARIZE_LEFT__HPP__ 1

#include <cicada/binarize_base.hpp>

#include <utils/indexed_set.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/dense_hash_map.hpp>

#include <boost/fusion/tuple.hpp>

namespace cicada
{
  struct BinarizeLeft : public BinarizeBase
  {
    typedef rule_type::symbol_set_type                 symbol_set_type;
    typedef hypergraph_type::edge_type::node_set_type  tail_set_type;

    template <typename Seq>
    struct hash_sequence : utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const Seq& x) const
      {
	return hasher_type::operator()(x.begin(), x.end(), 0);
      }
    };
    
    typedef utils::indexed_set<tail_set_type, hash_sequence<tail_set_type>, std::equal_to<tail_set_type>,
			       std::allocator<tail_set_type> > tail_map_type;
    typedef utils::indexed_set<symbol_set_type, hash_sequence<symbol_set_type>, std::equal_to<symbol_set_type>,
			       std::allocator<symbol_set_type> > symbol_map_type;
    
    typedef boost::fusion::tuple<tail_map_type::index_type, symbol_map_type::index_type, symbol_type> internal_type;
    typedef utils::dense_hash_map<internal_type, hypergraph_type::id_type, utils::hashmurmur<size_t>, std::equal_to<internal_type> >::type node_map_type;

    
    BinarizeLeft(const int __order=-1)
      : node_map(), order(__order)
    {
      node_map.set_empty_key(internal_type(-1, -1, symbol_type()));
    }
        
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

      if (! source.is_valid()) return;
      
      context_type context;
      phrase_type lhs_binarized;
      phrase_type binarized(2);
      hypergraph_type::edge_type::node_set_type tails(2);
      
      // we will traverse source-side in order to avoid confusion with newly created nodes...
      removed_type removed(source.edges.size(), false);

      position_set_type positions;
      
      tail_map.clear();
      symbol_map.clear();
      node_map.clear();

      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  
	  if (edge_source.tails.size() <= 2) continue;
	  
	  removed[edge_source.id] = true;
	  
	  // left most antecedents binarization...
	  
	  const symbol_type& lhs = edge_source.rule->lhs;
	  rule_type::symbol_set_type                rhs_sorted(edge_source.rule->rhs);
	  hypergraph_type::edge_type::node_set_type tails_sorted(edge_source.tails);
	  
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
	  
	  context.clear();
	  lhs_binarized.clear();
	  lhs_binarized.push_back(lhs);
	  for (size_t length = positions.size(); length > 2; -- length) {
	    const size_t first = 0;
	    const size_t last = first + length;
	    
	    context.push_back(rhs_sorted[positions[last - 1]].non_terminal());
	    
	    lhs_binarized.push_back(order < 0
				    ? binarized_label(lhs, context.begin(), context.end())
				    : binarized_label(lhs, std::max(context.begin(), context.end() - order), context.end()));
	  }
	  lhs_binarized.push_back(rhs_sorted[positions.front()]);
	  
	  // bottom-up traversal...
	  hypergraph_type::id_type tail_prev = tails_sorted.front();
	  phrase_type::const_reverse_iterator biter = lhs_binarized.rbegin();
	  for (size_t length = 2; length <= positions.size(); ++ length, ++ biter) {
	    const bool is_root = length == positions.size();
	    const size_t first = 0;
	    const size_t last = first + length;
	    const size_t middle = last - 1;

	    const symbol_type rhs_non_terminal = *biter;
	    const symbol_type lhs_non_terminal = *(biter + 1);
	    
	    binarized.clear();
	    
	    const size_t prefix_first = (is_root ? 0 : positions[first]);
	    const size_t prefix_last  = positions[first];
	    
	    binarized.insert(binarized.end(), rhs_sorted.begin() + prefix_first, rhs_sorted.begin() + prefix_last);
	    binarized.push_back(rhs_non_terminal);
	    
	    const size_t middle_first = positions[middle - 1] + 1;
	    const size_t middle_last  = positions[middle];
	    
	    binarized.insert(binarized.end(), rhs_sorted.begin() + middle_first, rhs_sorted.begin() + middle_last);
	    binarized.push_back(rhs_sorted[positions[last - 1]]);
	    
	    const size_t suffix_first = positions[last - 1] + 1;
	    const size_t suffix_last  = (is_root ? static_cast<int>(rhs_sorted.size()) : positions[last - 1] + 1);
	    
	    binarized.insert(binarized.end(), rhs_sorted.begin() + suffix_first, rhs_sorted.begin() + suffix_last);
	    
	    tails.front() = tail_prev;
	    tails.back()  = tails_sorted[last - 1];

	    if (is_root) {
	      hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	      edge_new.rule = rule_type::create(rule_type(lhs_non_terminal, binarized.begin(), binarized.end()));
	      
	      edge_new.features   = edge_source.features;
	      edge_new.attributes = edge_source.attributes;
	      
	      target.connect_edge(edge_new.id, node_source.id);
	    } else {
	      symbol_map_type::iterator siter = symbol_map.insert(symbol_set_type(binarized.begin(), binarized.end())).first;
	      tail_map_type::iterator   titer = tail_map.insert(tails).first;
	      
	      std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(internal_type(siter - symbol_map.begin(),
													     titer - tail_map.begin(),
													     lhs_non_terminal),
											       0));
	      
	      if (result.second) {
		const hypergraph_type::id_type node_id = target.add_node().id;
		
		hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
		
		edge_new.rule = rule_type::create(rule_type(lhs_non_terminal, binarized.begin(), binarized.end()));
		
		target.connect_edge(edge_new.id, node_id);
		
		result.first->second = edge_new.id;
	      }
	      
	      tail_prev = target.edges[result.first->second].head;
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
    
    tail_map_type   tail_map;
    symbol_map_type symbol_map;
    node_map_type   node_map;
    
    int order;
  };

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
