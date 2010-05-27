// -*- mode: c++ -*-

#ifndef __CICADA__INSIDE_OUTSIDE__HPP__
#define __CICADA__INSIDE_OUTSIDE__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>

#include <cicada/semiring/traits.hpp>

//
// inside-outside algorithms desribed in
//
//@InProceedings{li-eisner:2009:EMNLP,
//    author    = {Li, Zhifei  and  Eisner, Jason},
//    title     = {First- and Second-Order Expectation Semirings with Applications to Minimum-Risk Training on Translation Forests},
//    booktitle = {Proceedings of the 2009 Conference on Empirical Methods in Natural Language Processing},
//    month     = {August},
//    year      = {2009},
//    address   = {Singapore},
//    publisher = {Association for Computational Linguistics},
//    pages     = {40--51},
//    url       = {http://www.aclweb.org/anthology/D/D09/D09-1005}
//}
//

namespace cicada
{
  template <typename WeightSet, typename Function>
  inline
  void inside(const HyperGraph& graph, WeightSet& weights, Function function)
  {
    typedef typename WeightSet::value_type weight_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    weights.clear();
    weights.resize(graph.nodes.size(), weight_type());
    
    // visit in topological order... (we assume that the graph is "always" topologically ordered)
    hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
      const node_type& node = *niter;
      
      weight_type& weight = weights[node.id];
      
      node_type::edge_set_type::const_iterator eiter_end = node.in_edges.end();
      for (node_type::edge_set_type::const_iterator eiter = node.in_edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = graph.edges[*eiter];
	
	weight_type score = function(edge);
	edge_type::node_set_type::const_iterator niter_end = edge.tail_nodes.end();
	for (edge_type::node_set_type::const_iterator niter = edge.tail_nodes.begin(); niter != niter_end; ++ niter)
	  score *= weights[*niter];
	
	weight += score;
      }
    }
  };
  
  template <typename WeightSet, typename Function>
  inline
  void outside(const HyperGraph& graph, const WeightSet& weights_inside, WeightSet& weights_outside, Function function)
  {
    typedef typename WeightSet::value_type weight_type;
    
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    weights_outside.clear();
    weights_outside.resize(graph.nodes.size(), weight_type());

    weights_outside.back() = semiring::traits<weight_type>::one();

    if (weights_inside.size() != weights_outside.size())
      throw std::runtime_error("different inside/outside scores?");
    
    // visit in reversed topological order... (we assume that the graph is "always" topologically ordered)
    hypergraph_type::node_set_type::const_reverse_iterator niter_end = graph.nodes.rend();
    for (hypergraph_type::node_set_type::const_reverse_iterator niter = graph.nodes.rbegin(); niter != niter_end; ++ niter) {
      const node_type& node = *niter;

      const weight_type& score_head = weights_outside[node.id];
      
      node_type::edge_set_type::const_iterator eiter_end = node.in_edges.end();
      for (node_type::edge_set_type::const_iterator eiter = node.in_edges.begin(); eiter != eiter_end; ++ eiter) {
	
	const edge_type& edge = graph.edges[*eiter];
	
	weight_type score_head_edge = function(edge);
	score_head_edge *= score_head;
	
	edge_type::node_set_type::const_iterator niter_begin = edge.tail_nodes.begin();
	edge_type::node_set_type::const_iterator niter_end = edge.tail_nodes.end();
	for (edge_type::node_set_type::const_iterator niter = niter_begin; niter != niter_end; ++ niter) {
	  
	  weight_type score_outside = score_head_edge;
	  for (edge_type::node_set_type::const_iterator iiter = niter_begin; iiter != niter_end; ++ iiter)
	    if (iiter != niter)
	      score_outside *= weights_inside[*iiter];
	  
	  weights_outside[*niter] += score_outside;
	}
      }
    }
  }

  template <typename KWeight,   typename XWeight,
	    typename KFunction, typename XFunction>
  inline
  void inside_outside(const HyperGraph& graph,
		      KWeight& k,
		      XWeight& x,
		      KFunction function_k,
		      XFunction function_x)
  {
    typedef HyperGraph            hypergraph_type;

    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    typedef std::vector<KWeight, std::allocator<KWeight> > k_weight_set_type;

    k_weight_set_type inside_k;
    k_weight_set_type outside_k;
    
    inside<k_weight_set_type, KFunction>(graph, inside_k, function_k);
    outside<k_weight_set_type, KFunction>(graph, inside_k, outside_k, function_k);
    
    k = inside_k[graph.goal];
    x = XWeight();
    
    hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
      const node_type& node = *niter;
      
      node_type::edge_set_type::const_iterator eiter_end = node.in_edges.end();
      for (node_type::edge_set_type::const_iterator eiter = node.in_edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = graph.edges[*eiter];
	
	KWeight score_k = outside_k[node.id];
	
	edge_type::node_set_type::const_iterator niter_end = edge.tail_nodes.end();
	for (edge_type::node_set_type::const_iterator niter = edge.tail_nodes.begin(); niter != niter_end; ++ niter)
	  score_k *= inside_k[*niter];
	
	x += function_x(edge) * score_k;
      }
    }
  }
  
};

#endif
