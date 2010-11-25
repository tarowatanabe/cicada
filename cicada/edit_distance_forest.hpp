// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EDIT_DISTANCE_FOREST__HPP__
#define __CICADA__EDIT_DISTANCE_FOREST__HPP__ 1

#include <set>
#include <vector>

#include <cicada/symbol.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/vocab.hpp>

#include <utils/vector2.hpp>
#include <utils/bithack.hpp>

namespace cicada
{
  
  class EditDistanceForest
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;
    
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::id_type id_type;

  public:
    struct TRANSITION
    {
      enum transition_type {
	match,
	substitution,
	insertion,
	deletion,
      };
    };
    
  public:
    EditDistanceForest(const double __insertion = 1.0,
		       const double __deletion = 1.0,
		       const double __substitution = 1.0)
      : insertion(__insertion),
	deletion(__deletion),
	substitution(__substitution) {}
    
  private:
    typedef utils::vector2<double, std::allocator<double> > cost_matrix_type;
    
    typedef std::vector<id_type, std::allocator<id_type> >             node_set_type;
    
    typedef std::set<id_type, std::less<id_type>, std::allocator<id_type> > precedent_type;
    typedef std::vector<precedent_type, std::allocator<precedent_type> > lattice_type;
    
    typedef std::set<id_type, std::less<id_type>, std::allocator<id_type> > leaf_set_type;
    typedef std::vector<leaf_set_type, std::allocator<leaf_set_type> >      leaf_map_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > category_set_type;
    
  private:
    double insertion;
    double deletion;
    double subsitution;

    leaf_map_type leaf1;
    leaf_map_type leaf2;
    
    lattice_type lattice1;
    lattice_type lattice2;

    category_set_type cats1;
    category_set_type cats2;

    cost_matrix_type td;
    
  public:
    
    double operator()(const hypergraph_type& graph1, const hypergraph_tyep& graph2) const
    {
      if (! graph1.is_valid() && ! graph2.is_valid())
	return 0;
      else if (! graph1.is_valid())
	return graph2.nodes.size();
      else if (! graph2.is_valid())
	return graph1.nodes.size();

      td.clear();
      td.reserve(graph1.nodes.size(), graph2.nodes.size());
      td.resize(graph1.nodes.size(), graph2.nodes.size(), std::numeric_limits<double>::infinity());

      compute_category(graph1, cats1);
      compute_category(graph2, cats2);
      
      compute_left_most_leaf_descendant(graph1, leaf1);
      compute_left_most_leaf_descendant(graph2, leaf2);
      
      compute_precedent_nodes(graph1, leaf1, lattice1);
      compute_precedent_nodes(graph2, leaf2, lattice2);
      
      node_set_type key_roots1;
      node_set_type key_roots2;
      
      compute_key_roots(graph1, leaf1, key_roots1);
      compute_key_roots(graph2, leaf2, key_roots2);
      
      for (int x = 0; x != key_roots1.size(); ++ x)
	for (int y = 0; y != key_roots2.size(); ++ y)
	  for (leaf_set_type::const_iterator liter1 = leaf1[key_roots1[x]].begin(); liter1 != leaf1[key_roots1[x]].end(); ++ liter1)
	    for (leaf_set_type::const_iterator liter2 = leaf2[key_roots2[y]].begin(); liter2 != leaf2[key_roots2[y]].end(); ++ liter2)
	      forest_distance(graph1, graph2, *liter1, *liter2, key_roots1[x], key_roots2[y]);
      
      return td(graph1.goal, graph2.goal);
    }
    
  private:
    
    void forest_distance(const hypergraph_type& graph1,
			 const hypergraph_type& graph2,
			 const id_type& l1,
			 const id_type& l2,
			 const id_type& i,
			 const id_type& j)
    {
      cost_matrix_type fd(i + 2, j + 2, std::numeric_limits<double>::infinity());
      
      fd(l1, l2) = 0.0;
      
      for (id_type di = l1; di <= i; ++ di) {
	if (di == l1)
	  fd(di + 1, l2) = fd(di, l2) + deletion;
	else {
	  precedent_type::const_iterator piter_end = lattice1[di].end();
	  for (precedent_type::const_iterator piter = lattice1[di].begin(); piter != piter_end; ++ piter)
	    fd(di + 1, l2) = fd(*piter + 1, l2) + deletion;
	}
      }
      
      for (id_type dj = l2; dj <= j; ++ dj) {
	if (dj == l2)
	  fd(l1, dj + 1) = fd(l1, dj) + insertion;
	else {
	  precedent_type::const_iterator piter_end = lattice2[dj].end();
	  for (precedent_type::const_iterator piter = lattice2[dj].begin(); piter != piter_end; ++ piter)
	    fd(l1, dj + 1) = fd(l1, *piter + 1) + insertion;
	}
      }
      
      for (id_type di = l1; di <= i; ++ di)
	for (id_type dj = l2; dj <= j; ++ dj) {
	  if (leaf1[di].find(l1) != leaf1[di].end() && leaf2[dj].find(l2) != leaf2[dj].end()) {
	    // share the same l1/l2...

	    double& cost = fd(di + 1, dj + 1);
	    
	    // deletion...
	    if (di == l1)
	      cost = std::min(cost, fd(di, dj + 1) + deletion);
	    else {
	      precedent_type::const_iterator piter_end = lattice1[di].end();
	      for (precedent_type::const_iterator piter = lattice1[di].begin(); piter != piter_end; ++ piter)
		cost = std::min(cost, fd(*piter + 1, dj + 1) + deletion);
	    }

	    // insertion...
	    if (dj == l2)
	      cost = std::min(cost, fd(di + 1, dj) + insertion);
	    else {
	      precedent_type::const_iterator piter_end = lattice2[dj].end();
	      for (precedent_type::const_iterator piter = lattice2[dj].begin(); piter != piter_end; ++ piter)
		cost = std::min(cost, fd(di + 1, *piter + 1) + insertion);
	    }
	    
	    // rename...
	    const double cost_rename = (cats1[di] == cats2[dj] ? 0.0 : substitution);
	    
	    if (di == l1 && dj == l2)
	      cost = std::min(cost, fd(di, dj) + cost_rename);
	    else if (di == l1) {
	      precedent_type::const_iterator piter_end = lattice2[dj].end();
	      for (precedent_type::const_iterator piter = lattice2[dj].begin(); piter != piter_end; ++ piter)
		cost = std::min(cost, fd(di, *piter + 1) + cost_rename);
	    } else if (dj == l2) {
	      precedent_type::const_iterator piter_end = lattice1[di].end();
	      for (precedent_type::const_iterator piter = lattice1[di].begin(); piter != piter_end; ++ piter)
		cost = std::min(cost, fd(*piter + 1, dj) + cost_rename);
	    } else {
	      precedent_type::const_iterator piter1_end = lattice1[di].end();
	      precedent_type::const_iterator piter2_end = lattice2[dj].end();
	      for (precedent_type::const_iterator piter1 = lattice1[di].begin(); piter1 != piter1_end; ++ piter1)
		for (precedent_type::const_iterator piter2 = lattice2[dj].begin(); piter2 != piter2_end; ++ piter2)
		  cost = std::min(cost, fd(*piter1 + 1, *piter2 + 1) + cost_rename);
	    }
	    
	    td(di, dj) = std::min(td(di, dj), fd(di + 1, dj + 1));
	  } else {
	    // we do not share the same l1/l2...
	    
	    double& cost = fd(di + 1, dj + 1);
	    
	    // deletion...
	    if (di == l1)
	      cost = std::min(cost, fd(di, dj + 1) + deletion);
	    else {
	      precedent_type::const_iterator piter_end = lattice1[di].end();
	      for (precedent_type::const_iterator piter = lattice1[di].begin(); piter != piter_end; ++ piter)
		cost = std::min(cost, fd(*piter + 1, dj + 1) + deletion);
	    }
	    
	    // insertion...
	    if (dj == l2)
	      cost = std::min(cost, fd(di + 1, dj) + insertion);
	    else {
	      precedent_type::const_iterator piter_end = lattice2[dj].end();
	      for (precedent_type::const_iterator piter = lattice2[dj].begin(); piter != piter_end; ++ piter)
		cost = std::min(cost, fd(di + 1, *piter + 1) + insertion);
	    }
	    
	    
	    // enumerate leafs...
	    const double& cost_td = td(di, dj)
	    
	    leaf_set_type::const_iterator liter1_end = leaf1[di].end();
	    leaf_set_type::const_iterator liter2_end = leaf2[dj].end();
	    for (leaf_set_type::const_iterator liter1 = leaf1[di].begin(); liter1 != liter1_end; ++ liter1)
	      if (*liter1 >= l1)
		for (leaf_set_type::const_iterator liter2 = leaf2[dj].begin(); liter2 != liter2_end; ++ liter2)
		  if (*liter2 >= l2) {
		    
		    if (*liter1 == l1 && *liter2 == l2)
		      cost = std::min(cost, fd(*liter1, *liter2) + cost_td);
		    else if (*liter1 == l1) {
		      precedent_type::const_iterator piter_end = lattice2[*liter2].end();
		      for (precedent_type::const_iterator piter = lattice2[*liter2].begin(); piter != piter_end; ++ piter)
			cost = std::min(cost, fd(*liter1, *piter + 1) + cost_td);
		    } else if (*liter2 == l2) {
		      precedent_type::const_iterator piter_end = lattice1[*liter1].end();
		      for (precedent_type::const_iterator piter = lattice1[*liter1].begin(); piter != piter_end; ++ piter)
			cost = std::min(cost, fd(*piter + 1, *liter2) + cost_td);
		    } else {
		      precedent_type::const_iterator piter1_end = lattice1[*liter1].end();
		      precedent_type::const_iterator piter2_end = lattice2[*liter2].end();
		      for (precedent_type::const_iterator piter1 = lattice1[*liter1].begin(); piter1 != piter1_end; ++ piter1)
			for (precedent_type::const_iterator piter2 = lattice2[*liter2].begin(); piter2 != piter2_end; ++ piter2)
			  cost = std::min(cost, fd(*piter1 + 1, *piter2 + 1) + cost_td);
		    }
		  }
	    
	  }
	}
    }
    
    void compute_category(const hypergraph_type& graph, category_set_type& cats)
    {
      cats.clear();
      cats.reserve(graph.nodes.size());
      cats.resize(graph.nodes.size());
      
      for (id_type id = 0; id < graph.nodes.size(); ++ id) {
	const hypergraph_type::node_type& node = graph.nodes[id];
	
	if (node.edges.empty()) continue;
	
	const hypergraph_type::edge_type& edge = graph.edges[node.edges.front()];
	
	cats[id] = edge.rule->lhs;
      }
      
    }
    
    void compute_key_roots(const hypergraph_type& graph, const leaf_map_type& leaf, node_set_type& key_roots)
    {
      key_roots.clear();

      std::vector<bool, std::allocator<bool> > visited(graph.nodes.size(), false);
      
      for (int id = graph.goal; id >= 0; -- id) {
	
	bool inserted = false;
	
	leaf_set_type::const_iterator liter_end = leaf[id].end();
	for (leaf_set_type::const_iterator liter = leaf[id].begin(); liter != liter_end; ++ liter)
	  if (! visited[*liter]) {
	    if (! inserted)
	      key_roots.push_back(id);
	    inserted = true;
	    visited[*liter] = true;
	  }
      }
      
      std::sort(key_roots.begin(), key_roots.end());
    }
    
    void compute_left_most_leaf_descendant(const hypergraph_type& graph, leaf_map_type& leaf)
    {
      leaf.clear();
      leaf.reserve(graph.nodes.size());
      leaf.resize(graph.nodes.size());

      for (id_type id = 0; id != graph.nodes.size(); ++ id) {
	const hypergraph_type::node_type& node = graph.nodes[id];
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  if (! edge.tails.empty())
	    leaf[id].insert(leaf[edge.tails.front()].begin(), leaf[edge.tails.front()].end());
	  else
	    leaf[id].insert(id);
	}
      }
    }

    void compute_precedent_nodes(const hypergraph_type& graph, const leaf_map_type& leaf, lattice_type& lattice)
    {
      lattice.clear();
      lattice.reserve(graph.nodes.size());
      lattice.resize(graph.nodes.size());
      
      // topological order...? or, do we iterate from goal node?
      for (id_type id = 0; id != graph.nodes.size(); ++ id) {
	const hypergraph_type::node_type& node = graph.nodes[id];
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  if (! edge.tails.empty()) {
	    for (int i = 0; i < edge.tails.size() - 1; ++ i) {
	      // edge.tails[i] precedes nodes in leaf[edge.tails[i + 1]]
	      
	      leaf_set_type::const_iterator liter_end = leaf[edge.tails[i + 1]].end();
	      for (leaf_set_type::const_iterator liter = leaf[edge.tails[i + 1]].begin(); liter != liter_end; ++ liter)
		lattice[*liter].insert(edge.tails[i]);
	    }
	    
	    // edge.tails.back() precedes edge.head
	    lattice[edge.head].insert(edge.tails.back());
	  }
	}
      }
    }
  };
  
};

#endif
