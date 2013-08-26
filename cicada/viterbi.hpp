// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
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

    typedef Traversal traversal_type;
    typedef Function  function_type;
    
    typedef typename Traversal::value_type yield_type;
    typedef typename Function::value_type  semiring_type;
    typedef typename Function::value_type  weight_type;

    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
    typedef std::vector<yield_type, std::allocator<yield_type> >   derivation_set_type;
    
    typedef Viterbi<Traversal, Function> self_type;

    Viterbi(const hypergraph_type& __graph,
	    const traversal_type& __traversal,
	    const function_type& __function)
      : traversal(__traversal),
	function(__function),
	graph(__graph),
	weights(__graph.nodes.size()),
	derivations(__graph.nodes.size())
    {
      if (graph.goal == hypergraph_type::invalid)
	throw std::runtime_error("invalid hypergraph...");
    }
    
    friend struct Iterator;
    
    struct Iterator
    {
    private:
      typedef Viterbi<Traversal, Function> viterbi_type;

    public:
      typedef std::pair<weight_type, yield_type> value_type;
      typedef value_type* pointer;

    public:
      Iterator() : value(), viterbi(0) {}
      Iterator(const viterbi_type& __viterbi)
	: value(), viterbi(&const_cast<viterbi_type&>(__viterbi))
      {
	if (! viterbi->operator()(0, value.second, value.first)) {
	  value = value_type();
	  viterbi = 0;
	}
      }
      
    public:
      const value_type& operator*() const { return value; }
      const value_type* operator->() const { return &value; }
      
      Iterator& operator++()
      {
	if (viterbi) {
	  value = value_type();
	  viterbi = 0;
	}
	
	return *this;
      }

      Iterator operator++(int)
      {
	Iterator tmp = *this;
	++ *this;
	return tmp;
      }

      friend
      bool operator==(const Iterator& x, const Iterator& y)
      {
	return x.viterbi == y.viterbi;
      }

      friend
      bool operator!=(const Iterator& x, const Iterator& y)
      {
	return x.viterbi != y.viterbi;
      }

    private:
      value_type    value;
      viterbi_type* viterbi;
    };

    typedef Iterator iterator;
    typedef Iterator const_iterator;

  private:
    
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

  public:
    const_iterator begin() const { return const_iterator(*this); }
    iterator       begin() { return const_iterator(*this); }
    
    const_iterator end() const { return const_iterator(); }
    iterator       end() { return const_iterator(); }

  public:
    
    // k-best interface!
    bool operator()(int k, yield_type& yield, weight_type& weight)
    {
      // k is simply a dummy...
      
      yield_set_type yields;
      
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
      
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	
	  weight_type score = function(edge);
	  yields.clear();
	
	  // *=
	  hypergraph_type::edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	  for (hypergraph_type::edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter) {
	    score *= weights[*niter];
	    yields.push_back(&derivations[*niter]);
	  }
	  
	  // +=
	  if (score > weights[node.id]) {
	    weights[node.id] = score;
	    traversal(edge, derivations[node.id], yield_iterator(yields.begin()), yield_iterator(yields.end()));
	  }
	}
      }
      
      yield = derivations.back();
      weight = weights.back();

      return true;
    }

  private:
    const traversal_type traversal;
    const function_type  function;
    
    const hypergraph_type& graph;
    
    weight_set_type     weights;
    derivation_set_type derivations;
  };

  template <typename Traversal, typename Function>
  inline
  void viterbi(const HyperGraph& graph,
	       typename Traversal::value_type& yield,
	       typename Function::value_type&  weight,
	       const Traversal traversal,
	       const Function function)
  {
    Viterbi<Traversal, Function>(graph, traversal, function)(0, yield, weight);
  }
  
};

#endif
