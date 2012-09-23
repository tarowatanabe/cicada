// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_DEPENDENCY__HPP__
#define __CICADA__BINARIZE_DEPENDENCY__HPP__ 1

#include <cicada/binarize_base.hpp>

// we assume that every hyperedge takes the structure of:
// 
// [x]([y]^* terminal [z]^*)
//
// it is binarized so that:
// [x]([x^L] terminal [x^R])
// becomes
// [x^L]([x^L][y]) // left heavy...
// [x^L]([y] [y'])
// [x^R]([z][x^R]) // right heavy...
// [x^R]([z] [z'])

namespace cicada
{
  
  struct BinarizeDependency : public BinarizeBase
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
  };

  inline
  void binarize_dependency(const HyperGraph& source, HyperGraph& target)
  {
    BinarizeDependency binarizer;
    
    binarizer(source, target);
  }

  inline
  void binarize_dependency(HyperGraph& source)
  {
    HyperGraph target;

    BinarizeDependency binarizer;
    
    binarizer(source, target);
    
    source.swap(target);
  }

};


#endif
