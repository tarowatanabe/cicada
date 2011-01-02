// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__VITERBI__HPP__
#define __CICADA__VITERBI__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>

#include <cicada/semiring/traits.hpp>

namespace cicada
{
  // I tried semiring based implementation, but ended-up with specific function...

  template <typename Traversal, typename Function>
  struct Viterbi
  {
    typedef HyperGraph hypergraph_type;
    
    typedef typename Traversal::value_type yield_type;
    typedef typename Function::value_type  weight_type;
    
    typedef std::vector<const yield_type*, std::allocator<const yield_type*> > yield_set_type;
    
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
  };

  template <typename Traversal, typename Function>
  inline
  void viterbi(const HyperGraph& graph,
	       typename Traversal::value_type& yield,
	       typename Function::value_type&  weight,
	       const Traversal traversal,
	       const Function function)
  {
    typedef HyperGraph hypergraph_type;

    typedef typename Traversal::value_type yield_type;
    typedef typename Function::value_type  weight_type;
    
    typedef std::vector<yield_type, std::allocator<yield_type> >   yield_set_type;
    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;

    typedef std::vector<const yield_type*, std::allocator<const yield_type*> > yield_antecedent_type;

    yield_set_type yields(graph.nodes.size());
    weight_set_type weights(graph.nodes.size());

    yield_antecedent_type antecedents;
    
    hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
      const hypergraph_type::node_type& node = *niter;
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	
	weight_type score = function(edge);
	
	antecedents.clear();
	
	// *=
	hypergraph_type::edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	for (hypergraph_type::edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter) {
	  score *= weights[*niter];
	  antecedents.push_back(&yields[*niter]);
	}

	typedef typename Viterbi<Traversal, Function>::yield_iterator yield_iterator;
	
	// +=
	if (score > weights[node.id]) {
	  weights[node.id] = score;
	  traversal(edge, yields[node.id], yield_iterator(antecedents.begin()), yield_iterator(antecedents.end()));
	}
      }
    }
    
    yield = yields.back();
    weight = weights.back();
  }
  
};

#endif
