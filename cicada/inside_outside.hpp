// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

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
  
  template <typename Function>
  struct Inside
  {
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    Inside(Function __function) : function(__function) {}

    template <typename WeightSet>
    void operator()(const hypergraph_type& graph, WeightSet& weights)
    {
      typedef typename WeightSet::value_type weight_type;
      
      // visit in topological order... (we assume that the graph is "always" topologically ordered)
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
	const node_type& node = *niter;
	
	weight_type& weight = weights[node.id];
	
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const edge_type& edge = graph.edges[*eiter];
	  
	  weight_type score = function(edge);
	  edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	  for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	    score *= weights[*niter];
	  
	  weight += score;
	}
      }
    }
    
    Function function;
  };

  template <typename Function>
  struct Outside
  {
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    Outside(Function __function) : function(__function) {}
    
    template <typename WeightSet, typename WeightSetOutside>
    void operator()(const hypergraph_type& graph, const WeightSet& weights_inside, WeightSetOutside& weights_outside)
    {
      typedef typename WeightSet::value_type weight_type;
      
      weights_outside[graph.nodes.size() - 1] = semiring::traits<weight_type>::one();

      if (weights_inside.size() != weights_outside.size())
	throw std::runtime_error("different inside/outside scores?");
    
      // visit in reversed topological order... (we assume that the graph is "always" topologically ordered)
      hypergraph_type::node_set_type::const_reverse_iterator niter_end = graph.nodes.rend();
      for (hypergraph_type::node_set_type::const_reverse_iterator niter = graph.nodes.rbegin(); niter != niter_end; ++ niter) {
	const node_type& node = *niter;

	const weight_type& score_head = weights_outside[node.id];
      
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	
	  const edge_type& edge = graph.edges[*eiter];
	
	  weight_type score_head_edge = function(edge);
	  score_head_edge *= score_head;
	
	  edge_type::node_set_type::const_iterator niter_begin = edge.tails.begin();
	  edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
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
    
    Function function;
  };
  
  template <typename KFunction, typename XFunction>
  struct InsideOutside
  {
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    InsideOutside(KFunction __function_k,
		  XFunction __function_x)
      : inside(__function_k),
	outside(__function_k),
	function_x(__function_x) {}

    template <typename KWeightSet, typename KWeightOutsideSet, typename XWeightSet>
    void operator()(const hypergraph_type& graph,
		    KWeightSet& inside_k,
		    KWeightOutsideSet& outside_k,
		    XWeightSet& x)
    {
      typedef typename KWeightSet::value_type KWeight;
      typedef typename XWeightSet::value_type XWeight;
      
      inside(graph, inside_k);
      outside(graph, inside_k, outside_k);
    
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
	const node_type& node = *niter;
      
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const edge_type& edge = graph.edges[*eiter];
	
	  KWeight score_k = outside_k[node.id];
	
	  edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	  for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	    score_k *= inside_k[*niter];
	  
	  x[edge.id] += function_x(edge) * score_k;
	}
      }
    }
    
    
    template <typename KWeightSet, typename XWeightSet>
    void operator()(const hypergraph_type& graph,
		    KWeightSet& inside_k,
		    XWeightSet& x)
    {
      typedef typename KWeightSet::value_type KWeight;
      typedef typename XWeightSet::value_type XWeight;
      
      typedef std::vector<KWeight, std::allocator<KWeight> > k_weight_set_type;
      
      k_weight_set_type outside_k(graph.nodes.size());
      
      operator()(graph, inside_k, outside_k, x);
    }
    
    Inside<KFunction>  inside;
    Outside<KFunction> outside;
    XFunction function_x;
  };

  template <typename WeightSet, typename Function>
  inline
  void inside(const HyperGraph& graph, WeightSet& weights, Function function)
  {
    Inside<Function> __inside(function);
    
    __inside(graph, weights);
  };
  
  template <typename WeightSet, typename WeightSetOutside, typename Function>
  inline
  void outside(const HyperGraph& graph, const WeightSet& weights_inside, WeightSetOutside& weights_outside, Function function)
  {
    Outside<Function> __outside(function);
    
    __outside(graph, weights_inside, weights_outside);
  }

  template <typename KWeightSet, typename XWeightSet,
	    typename KFunction, typename XFunction>
  inline
  void inside_outside(const HyperGraph& graph,
		      KWeightSet& inside_k,
		      XWeightSet& x,
		      KFunction function_k,
		      XFunction function_x)
  {
    InsideOutside<KFunction, XFunction> __inside_outside(function_k, function_x);
    
    __inside_outside(graph, inside_k, x);
  }
  
  template <typename KWeightSet, typename KWeightOutsideSet, typename XWeightSet,
	    typename KFunction, typename XFunction>
  inline
  void inside_outside(const HyperGraph& graph,
		      KWeightSet& inside_k,
		      KWeightOutsideSet& outside_k,
		      XWeightSet& x,
		      KFunction function_k,
		      XFunction function_x)
  {
    InsideOutside<KFunction, XFunction> __inside_outside(function_k, function_x);
    
    __inside_outside(graph, inside_k, outside_k, x);
  }
  
};

#endif
