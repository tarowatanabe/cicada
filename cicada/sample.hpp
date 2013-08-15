// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//
// a sampler based on the KBest-like interface

#ifndef __CICADA__SAMPLE__HPP__
#define __CICADA__SAMPLE__HPP__ 1

#include <vector>
#include <algorithm>

#include <cicada/hypergraph.hpp>
#include <cicada/inside_outside.hpp>
#include <cicada/semiring/traits.hpp>

#include <utils/bithack.hpp>

namespace cicada
{
  template <typename Traversal, typename Function, typename Sampler>
  struct Sample
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    typedef Traversal traversal_type;
    typedef Function  function_type;
    typedef Sampler   sampler_type;

    typedef typename traversal_type::value_type yield_type;
    typedef typename function_type::value_type  semiring_type;
    typedef typename function_type::value_type  weight_type;

    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;

    typedef std::vector<yield_type, std::allocator<yield_type> > derivation_set_type;
    
    Sample(const hypergraph_type& __graph,
	   const traversal_type& __traversal,
	   const function_type& __function,
	   sampler_type& __sampler,
	   const double __temperature=1.0)
      : traversal(__traversal),
	function(__function),
	sampler(__sampler),
	graph(__graph),
	insides(__graph.nodes.size()),
	derivations(__graph.nodes.size()),
	temperature(__temperature)
    {
      if (graph.goal == hypergraph_type::invalid)
	throw std::runtime_error("invalid hypergraph...");
      
      cicada::inside(graph, insides, function);
    }

    
    bool operator()(int k, yield_type& yield, weight_type& weight)
    {
      //
      // we will simply ignore k, but leave it as is to preserves the KBest interface
      //
      
      return sample_derivation(yield, weight);
    }

  private:
    typedef std::vector<const yield_type*, std::allocator<const yield_type*> > yield_set_type;

  public:
    class yield_iterator : public yield_set_type::const_iterator
    {
    public:
      typedef typename yield_set_type::const_iterator base_type;
      
      yield_iterator(const base_type& x) : base_type(x) {}
      
      const yield_type& operator*()  { return *(base_type::operator*()); }
      const yield_type* operator->() { return base_type::operator*(); }
      
      friend
      yield_iterator operator+(const yield_iterator& x, ptrdiff_t diff)
      {
	return yield_iterator(base_type(x) + diff);
      }

      friend
      yield_iterator operator-(const yield_iterator& x, ptrdiff_t diff)
      {
	return yield_iterator(base_type(x) - diff);
      }
    };

  private:
    typedef std::vector<id_type, std::allocator<id_type> > stack_type;
    typedef std::pair<id_type, id_type> node_edge_type;
    typedef std::vector<node_edge_type, std::allocator<node_edge_type> > node_edge_set_type;

    struct compare_node_edge
    {
      bool operator()(const node_edge_type& x, const node_edge_type& y) const
      {
	return x.first < y.first;
      }
    };

    stack_type         stack;
    node_edge_set_type edges;
    weight_set_type    scores;
    weight_set_type    probs;
    
    bool sample_derivation(yield_type& yield, weight_type& weight)
    {
      // perform top-down traversal for samling a tree

      weight = cicada::semiring::traits<weight_type>::one();
      
      stack.clear();
      stack.push_back(graph.goal);

      edges.clear();
      
      while (! stack.empty()) {
	const id_type node_id =  stack.back();
	stack.pop_back();
	
	const node_type& node = graph.nodes[node_id];

	// this is an error!
	if (node.edges.empty()) return false;
	
	size_type pos_sampled = 0;
	if (node.edges.size() == 1)
	  weight *= function(graph.edges[node.edges.front()]);
	else {
	  weight_type sum;
	  
	  scores.clear();
	  probs.clear();
	  
	  node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	  for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	    const edge_type& edge = graph.edges[*eiter];
	    
	    weight_type prob = function(edge);
	    
	    scores.push_back(prob);
	    
	    edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	    for (edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	      prob *= insides[*titer];
	    
	    sum += prob;
	    probs.push_back(prob);
	  }
	  
	  // normalize... if summation is zero, then, use uniform distribution!
	  if (sum != cicada::semiring::traits<weight_type>::zero())
	    std::transform(probs.begin(), probs.end(), probs.begin(), std::bind2nd(std::multiplies<weight_type>(), 1.0 / sum));
	  else
	    std::fill(probs.begin(), probs.end(), weight_type(1.0 / probs.size()));
	  
	  pos_sampled = sampler.draw(probs.begin(), probs.end(), temperature) - probs.begin();
	  
	  // updated weight...
	  weight *= scores[pos_sampled];
	}
	
	// sampled edge-id
	const id_type edge_id_sampled = node.edges[pos_sampled];
	
	// update stack...
	const edge_type& edge_sampled = graph.edges[edge_id_sampled];
	edge_type::node_set_type::const_iterator titer_end = edge_sampled.tails.end();
	for (edge_type::node_set_type::const_iterator titer = edge_sampled.tails.begin(); titer != titer_end; ++ titer)
	  stack.push_back(*titer);
	
	// update sampled edges
	edges.push_back(std::make_pair(node_id, edge_id_sampled));
      }

      if (edges.empty()) return false;
      
      // perform bottom-up to collect yields
      // First, we need to make sure that the edges will be visited in a topilogical order...
      std::sort(edges.begin(), edges.end(), compare_node_edge());
      
      // Second, collect yields via traversals
      yield_set_type yields;
      
      typename node_edge_set_type::const_iterator eiter_end = edges.end();
      for (typename node_edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = graph.edges[eiter->second];

	yields.clear();
	edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	for (edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	  yields.push_back(&derivations[*titer]);
	
	traversal(edge, derivations[eiter->first], yield_iterator(yields.begin()), yield_iterator(yields.end()));
      }
      
      // final yield!
      yield = derivations[graph.goal];

      return true;
    }


  private:
    const traversal_type traversal;
    const function_type  function;
    sampler_type&        sampler;
    
    const hypergraph_type& graph;
    
    weight_set_type     insides;
    derivation_set_type derivations;
    
    const double temperature;
  };
};

#endif
