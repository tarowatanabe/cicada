// -*- mode: c++ -*-
//
//  Copyright(C) 2012-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_DEPENDENCY__HPP__
#define __CICADA__BINARIZE_DEPENDENCY__HPP__ 1

#include <cicada/binarize_base.hpp>

#include <utils/hashmurmur3.hpp>

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

    typedef std::vector<tails_type, std::allocator<tails_type> > tails_map_type;
    typedef std::vector<rhs_type, std::allocator<rhs_type> > rhs_map_type;

    typedef rule_type::symbol_set_type                 symbol_set_type;
    typedef hypergraph_type::edge_type::node_set_type  tail_set_type;

    typedef std::pair<symbol_type, tail_set_type> binarized_type;
    
    struct binarized_hash : public utils::hashmurmur3<size_t>
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      size_t operator()(const binarized_type& x) const {
	return hasher_type::operator()(x.second.begin(), x.second.end(), x.first.id());
      }
    };
    
    typedef utils::unordered_map<binarized_type, hypergraph_type::id_type, binarized_hash, std::equal_to<binarized_type>,
				 std::allocator<std::pair<const binarized_type, hypergraph_type::id_type> > >::type binarized_map_type;

    BinarizeDependency(const bool __head_mode=false,
		       const bool __label_mode=false)
      : head_mode(__head_mode),
	label_mode(__label_mode) {}
    
    const bool head_mode;
    const bool label_mode;
    
    binarized_map_type binarized;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // first, copy...
      target = source;
      
      if (! source.is_valid()) return;
      
      rhs_map_type   rhs;
      tails_map_type tails;
      
      rhs_type   rhs_new;
      tails_type tails_new;
      
      removed_type removed(source.edges.size(), false);
      
      binarized.clear();
      
      const std::string tag_neutral("^");
      const std::string tag_left("^L");
      const std::string tag_right("^R");
      const std::string tag_middle("^M");
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  
	  // for head_mode mode, even single non-terminal should be checked!
	  if (edge_source.tails.size() < 2 - head_mode) continue;
	  
	  rhs.clear();
	  tails.clear();
	  
	  {
	    enum {
	      NONE = 0,
	      TERMINAL = 1,
	      NON_TERMINAL = 2,
	    };
	    
	    int mode = NONE;
	    int non_terminal_pos = 0;
	    
	    rule_type::symbol_set_type::const_iterator riter_end = edge_source.rule->rhs.end();
	    for (rule_type::symbol_set_type::const_iterator riter = edge_source.rule->rhs.begin(); riter != riter_end; ++ riter)
	      if (riter->is_non_terminal()) {
		const int __non_terminal_index = riter->non_terminal_index();
		const int non_terminal_index = utils::bithack::branch(__non_terminal_index <= 0,
								      non_terminal_pos,
								      __non_terminal_index - 1);
		++ non_terminal_pos;
		
		if (mode == NONE || mode == TERMINAL) {
		  rhs.push_back(rhs_type());
		  tails.push_back(tails_type());
		}
		mode = NON_TERMINAL;
		
		rhs.back().push_back(riter->non_terminal());
		tails.back().push_back(edge_source.tails[non_terminal_index]);
	      } else {
		if (mode == NONE || mode == NON_TERMINAL) {
		  rhs.push_back(rhs_type());
		  tails.push_back(tails_type());
		}
		mode = TERMINAL;
		
		rhs.back().push_back(*riter);
	      }
	  }
	  
	  const symbol_type& lhs = edge_source.rule->lhs;
	  
	  rhs_new.clear();
	  tails_new.clear();
	  
	  tails_map_type::const_iterator titer = tails.begin();
	  rhs_map_type::const_iterator riter_begin = rhs.begin();
	  rhs_map_type::const_iterator riter_end   = rhs.end();
	  for (rhs_map_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter, ++ titer) {
	    if (titer->empty()) {
	      if (head_mode) {
		const std::pair<symbol_type, hypergraph_type::id_type> result = headlize(lhs, *riter, target);
		
		rhs_new.push_back(result.first);
		tails_new.push_back(result.second);
	      } else
		rhs_new.insert(rhs_new.end(), riter->begin(), riter->end());
	    } else if (titer->size() == 1) {
	      rhs_new.push_back(riter->front());
	      tails_new.push_back(titer->front());
	    } else {
	      const std::string& tag = (! label_mode || rhs.size() == 1
					? tag_neutral
					: (riter == riter_begin
					   ? tag_left
					   : (riter + 1 == riter_end
					      ? tag_right
					      : tag_middle)));
	      
	      const std::pair<symbol_type, hypergraph_type::id_type> result = binarize(tag, *riter, *titer, target);
	      
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

    std::pair<symbol_type, hypergraph_type::id_type> headlize(const symbol_type& root,
							      const rhs_type& rhs,
							      hypergraph_type& target)
    {
      const symbol_type lhs = (label_mode
			       ? '[' + root.non_terminal_strip() + "^H]"
			       : '[' + root.non_terminal_strip() + "^]");
      
      const hypergraph_type::id_type head = target.add_node().id;
      
      hypergraph_type::edge_type& edge_new = target.add_edge();
      edge_new.rule = rule_type::create(rule_type(lhs, rhs.begin(), rhs.end()));
      target.connect_edge(edge_new.id, head);
      
      return std::make_pair(lhs, head);
    }
    
    std::pair<symbol_type, hypergraph_type::id_type> binarize(const std::string& tag,
							      const rhs_type& rhs,
							      const tails_type& tails,
							      hypergraph_type& target)
    {
      typedef size_t size_type;
      
      const size_type child_size = rhs.size();

      label_chart.clear();
      label_chart.resize(child_size + 1, symbol_type());
      
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

	  {
	    const symbol_type::piece_type left = label_chart(first, last - 1).non_terminal_strip();
	    const symbol_type::piece_type right = label_chart(last - 1, last).non_terminal_strip();
	    
	    if (length > 2)
	      label_chart(first, last) = '[' + std::string(left.begin(), left.end() - tag.size()) + '+' + right + tag + ']';
	    else
	      label_chart(first, last) = '[' + std::string(left) + '+' + right + tag + ']';
	  }
	  
	  const symbol_type lhs = label_chart(first, last);
	  
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
	    } else {
	      if (head == hypergraph_type::invalid)
		head = target.edges[result.first->second].head;
	      else if (head != target.edges[result.first->second].head)
		throw std::runtime_error("invalid binarized edge");
	    }
	  }
	}
      
      return std::make_pair(label_chart(0, child_size), node_chart(0, child_size));
    }
    
  };

  inline
  void binarize_dependency(const HyperGraph& source, HyperGraph& target, const bool head_mode=false, const bool label_mode=false)
  {
    BinarizeDependency binarizer(head_mode, label_mode);
    
    binarizer(source, target);
  }

  inline
  void binarize_dependency(HyperGraph& source, const bool head_mode=false, const bool label_mode=true)
  {
    HyperGraph target;

    BinarizeDependency binarizer(head_mode, label_mode);
    
    binarizer(source, target);
    
    source.swap(target);
  }

};


#endif
