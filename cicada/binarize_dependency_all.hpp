// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_DEPENDENCY_ALL__HPP__
#define __CICADA__BINARIZE_DEPENDENCY_ALL__HPP__ 1

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
  
  struct BinarizeDependencyAll : public BinarizeBase
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

    binarized_map_type binarized;
    
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

      binarized.clear();
      
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
	      const symbol_type lhs_left = "[" + lhs.non_terminal_strip() + "^L]";
	      
	      const std::pair<symbol_type, hypergraph_type::id_type> result = binarize(lhs_left, rhs_left, tails_left, target);
	      
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
	      const symbol_type lhs_right = "[" + lhs.non_terminal_strip() + "^R]";

	      const std::pair<symbol_type, hypergraph_type::id_type> result = binarize(lhs_right, rhs_right, tails_right, target);
	      
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

  private:
    node_chart_type   node_chart;
    label_chart_type  label_chart;
    
    std::pair<symbol_type, hypergraph_type::id_type> binarize(const symbol_type& lhs,
							      const rhs_type& rhs,
							      const tails_type& tails,
							      hypergraph_type& target)
    {
      typedef size_t size_type;
      
      const size_type child_size = rhs.size();

      label_chart.clear();
      label_chart.resize(child_size + 1, lhs);
      
      node_chart.clear();
      node_chart.resize(child_size + 1, hypergraph_type::invalid);
      
      for (size_type i = 0; i != child_size; ++ i) {
	label_chart(i, i + 1) = rhs[i];
	node_chart(i, i + 1)  = tails[i];
      }
      
      symbol_set_type rhs_binarized(2);
      tail_set_type   tails_binarized(2);
  
      for (size_type length = 2; length <= child_size; ++ length)
	for (size_type first = 0; first + length <= child_size; ++ first) {
	  const size_type last = first + length;
	  
	  hypergraph_type::id_type& head = node_chart(first, last);
	  
	  for (size_t middle = first + 1; middle != last; ++ middle) {
	    // [first, middle) and [middle, last)
	    
	    rhs_binarized.front() = label_chart(first, middle);
	    rhs_binarized.back()  = label_chart(middle, last);
	    
	    tails_binarized.front() = node_chart(first, middle);
	    tails_binarized.back()  = node_chart(middle, last);
	    
	    std::pair<binarized_map_type::iterator, bool> result = binarized.insert(std::make_pair(binarized_type(lhs,
														  tails_binarized),
												   0));
	    if (result.second) {
	      if (head == hypergraph_type::invalid)
		head = target.add_node().id;
	      
	      hypergraph_type::edge_type& edge_new = target.add_edge(tails_binarized.begin(), tails_binarized.end());
	      
	      edge_new.rule = rule_type::create(rule_type(lhs, rhs_binarized.begin(), rhs_binarized.end()));
	      
	      target.connect_edge(edge_new.id, head);
	      
	      result.first->second = edge_new.id;
	    } else
	      head = target.edges[result.first->second].head;
	  }
	}
      
      return std::make_pair(label_chart(0, child_size), node_chart(0, child_size));
    }
    
  };

  inline
  void binarize_dependency_all(const HyperGraph& source, HyperGraph& target)
  {
    BinarizeDependencyAll binarizer;
    
    binarizer(source, target);
  }

  inline
  void binarize_dependency_all(HyperGraph& source)
  {
    HyperGraph target;

    BinarizeDependencyAll binarizer;
    
    binarizer(source, target);
    
    source.swap(target);
  }

};


#endif
