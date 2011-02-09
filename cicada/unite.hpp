// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__UNITE__HPP__
#define __CICADA__UNITE__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/lattice.hpp>
#include <cicada/sort.hpp>

namespace cicada
{
  struct Unite
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    typedef Lattice lattice_type;

    typedef lattice_type::feature_set_type feature_set_type;

    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;
    typedef Rule   rule_type;

    
    void operator()(lattice_type& x, const lattice_type& y)
    {
      if (&x == &y || y.empty()) return;
      if (x.empty()) {
	x = y;
	return;
      }
      
      lattice_type merged;
      
      merged.push_back(lattice_type::arc_set_type());
      merged.back().push_back(lattice_type::arc_type(vocab_type::EPSILON, feature_set_type(), 1));
      merged.back().push_back(lattice_type::arc_type(vocab_type::EPSILON, feature_set_type(), x.size() + 2));
      
      for (size_t i = 0; i != x.size(); ++ i)
	merged.push_back(x[i]);
      
      merged.push_back(lattice_type::arc_set_type());
      merged.back().push_back(lattice_type::arc_type(vocab_type::EPSILON, feature_set_type(), y.size() + 1));
      
      for (size_t i = 0; i != y.size(); ++ i)
	merged.push_back(y[i]);
      
      x.swap(merged);

      x.initialize_distance();
    }
    
    void operator()(hypergraph_type& x, const hypergraph_type& y)
    {
      if (&x == &y || y.nodes.empty() || y.goal == hypergraph_type::invalid) return;
      if (x.nodes.empty() || x.goal == hypergraph_type::invalid) {
	x = y;
	return;
      }
      
      //
      // check if we share the same goal... otherwise, we will create new edge and fill...
      // 
      const symbol_type goal_x = x.edges[x.nodes[x.goal].edges.front()].rule->lhs;
      const symbol_type goal_y = y.edges[y.nodes[y.goal].edges.front()].rule->lhs;

      if (goal_x == goal_y) {
	const int y_node_offset = x.nodes.size();
	const int y_edge_offset = x.edges.size();
	
	// -1 to adjust goal
	x.nodes.resize(x.nodes.size() + y.nodes.size() - 1);
	x.edges.resize(x.edges.size() + y.edges.size());
	
	// merge nodes
	for (id_type id = 0; id < y.nodes.size(); ++ id)
	  if (id != y.goal) {
	    const id_type id_new = id + y_node_offset - (id > y.goal);
	    
	    const node_type& node_old = y.nodes[id];
	    node_type& node_new = x.nodes[id_new];
	    
	    node_new = node_old;
	    
	    node_new.id =  id_new;
	    node_type::edge_set_type::iterator eiter_end = node_new.edges.end();
	    for (node_type::edge_set_type::iterator eiter = node_new.edges.begin(); eiter != eiter_end; ++ eiter)
	      *eiter += y_edge_offset;
	  }
	
	// merge edges
	for (id_type id = 0; id < y.edges.size(); ++ id) {
	  const edge_type& edge_old = y.edges[id];
	  edge_type& edge_new = x.edges[id + y_edge_offset];
	  
	  edge_new = edge_old;
	  edge_new.id = id + y_edge_offset;
	  
	  if (edge_new.head == y.goal) {
	    edge_new.head = x.goal;
	    x.nodes[x.goal].edges.push_back(edge_new.id);
	  } else
	    edge_new.head += y_node_offset - (edge_new.head > y.goal);
	  
	  edge_type::node_set_type::iterator niter_end = edge_new.tails.end();
	  for (edge_type::node_set_type::iterator niter = edge_new.tails.begin(); niter != niter_end; ++ niter)
	    *niter += y_node_offset - (*niter > y.goal);
	}
      } else {
	const int y_node_offset = x.nodes.size();
	const int y_edge_offset = x.edges.size();
	
	// +1 to adjust goal
	x.nodes.resize(x.nodes.size() + y.nodes.size() + 1);
	// +2 to adjust edges toward goal
	x.edges.resize(x.edges.size() + y.edges.size() + 2);
	
	// merge nodes
	for (id_type id = 0; id < y.nodes.size(); ++ id) {
	  const node_type& node_old = y.nodes[id];
	  node_type& node_new = x.nodes[id + y_node_offset];
	  
	  node_new = node_old;
	    
	  node_new.id = id + y_node_offset;
	  node_type::edge_set_type::iterator eiter_end = node_new.edges.end();
	  for (node_type::edge_set_type::iterator eiter = node_new.edges.begin(); eiter != eiter_end; ++ eiter)
	    *eiter += y_edge_offset;
	}
	
	// merge edges
	for (id_type id = 0; id < y.edges.size(); ++ id) {
	  const edge_type& edge_old = y.edges[id];
	  edge_type& edge_new = x.edges[id + y_edge_offset];
	  
	  edge_new = edge_old;
	  edge_new.id = id + y_edge_offset;
	  
	  edge_new.head += y_node_offset;
	  edge_type::node_set_type::iterator niter_end = edge_new.tails.end();
	  for (edge_type::node_set_type::iterator niter = edge_new.tails.begin(); niter != niter_end; ++ niter)
	    *niter += y_node_offset;
	}
	
	// create new node...
	node_type& goal_node = x.nodes.back();
	edge_type& goal_edge_x = x.edges[x.edges.size() - 2];
	edge_type& goal_edge_y = x.edges[x.edges.size() - 1];
	
	goal_node.id = x.nodes.size() - 1;
	goal_edge_x.id = x.edges.size() - 2;
	goal_edge_y.id = x.edges.size() - 1;

	goal_edge_x.rule = rule_type::create(rule_type(vocab_type::GOAL,
						       rule_type::symbol_set_type(1, goal_x.non_terminal())));
	goal_edge_y.rule = rule_type::create(rule_type(vocab_type::GOAL,
						       rule_type::symbol_set_type(1, goal_y.non_terminal())));
	
	goal_edge_x.head = goal_node.id;
	goal_edge_y.head = goal_node.id;
	
	goal_edge_x.tails = edge_type::node_set_type(1, x.goal);
	goal_edge_y.tails = edge_type::node_set_type(1, y.goal + y_node_offset);
	
	goal_node.edges.push_back(goal_edge_x.id);
	goal_node.edges.push_back(goal_edge_y.id);
	
	x.goal = goal_node.id;
      }
      
      cicada::topologically_sort(x);
    }
  };

  inline
  void unite(const HyperGraph& x, const HyperGraph& y, HyperGraph& merged)
  {
    Unite __unite;
    
    merged = x;
    
    __unite(merged, y);
  }
  
  inline
  void unite(HyperGraph& x, const HyperGraph& y)
  {
    Unite __unite;
    
    __unite(x, y);
  }

  inline
  void unite(const Lattice& x, const Lattice& y, Lattice& merged)
  {
    Unite __unite;
    
    merged = x;
    
    __unite(merged, y);
  }
  
  inline
  void unite(Lattice& x, const Lattice& y)
  {
    Unite __unite;
    
    __unite(x, y);
  }
  
};



#endif
