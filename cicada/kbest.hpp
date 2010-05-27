// -*- mode: c++ -*-

#ifndef __CICADA__KBEST__HPP__
#define __CICADA__KBEST__HPP__ 1

#include <vector>
#include <algorithm>

#include <cicada/hypergraph.hpp>
#include <cicada/semiring/traits.hpp>

#include <utils/bithack.hpp>
#include <utils/simple_vector.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/hashmurmur.hpp>

namespace cicada
{
  // k-best derication based on Algorithm 3 of
  //
  // @InProceedings{huang-chiang:2005:IWPT,
  //  author    = {Huang, Liang  and  Chiang, David},
  //  title     = {Better k-best Parsing},
  //  booktitle = {Proceedings of the Ninth International Workshop on Parsing Technology},
  //  month     = {October},
  //  year      = {2005},
  //  address   = {Vancouver, British Columbia},
  //  publisher = {Association for Computational Linguistics},
  //  pages     = {53--64},
  //  url       = {http://www.aclweb.org/anthology/W/W05/W05-1506}
  //  }

  
  // semiring,
  // yield, requires operator=(const yield&) (assignment) and  yield::yield() (constructor)
  // traversal function operator()(const HyperGraph::Edge&, const yield*, Iterator first, Iterator last);
  //                    where Iterator's value (*first etc.) is const yield*
  // semiring function
  
  
  template <typename Semiring,
	    typename Yield,
	    typename Traversal,
	    typename Function>
  struct KBest
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    typedef Semiring semiring_type;
    typedef Semiring weight_type;
    typedef Yield    yield_type;
    typedef Traversal traversal_type;
    typedef Function  function_type;
    
    KBest(const hypergraph_type& __graph,
	  const size_type& __k_prime,
	  const traversal_type& __traversal,
	  const function_type& __function)
      : graph(__graph),
	k_prime(__k_prime),
	traversal(__traversal),
	function(__function),
	states(__graph.nodes.size()) {}
    
    typedef utils::simple_vector<int, std::allocator<int> > index_set_type;
    
    struct Derivation
    {
      Derivation(const index_set_type& __j) : j(__j) {}
      
      yield_type yield;
      const edge_type* edge;
      const index_set_type j;
      weight_type      score;
    };
    
    typedef Derivation derivation_type;
    typedef utils::chunk_vector<derivation_type, 4096 / sizeof(derivation_type), std::allocator<derivation_type> > derivation_set_type;

    struct compare_derivation_type
    {
      bool operator()(const derivation_type* x, const derivation_type* y) const
      {
	return x->score > y->score;
      }
    };
    
    struct compare_heap_type
    {
      bool operator()(const derivation_type* x, const derivation_type* y) const
      {
	return x->score < y->score;
      }
    };
    
    typedef std::vector<const derivation_type*, std::allocator<const derivation_type*> > derivation_heap_type;
    typedef std::vector<const derivation_type*, std::allocator<const derivation_type*> > derivation_list_type;
    
    struct derivation_hash_type : public utils::hashmurmur<size_t>
    {
      size_t operator()(const derivation_type* x) const
      {
	return utils::hashmurmur<size_t>::operator()(x->j.begin(), x->j.end(), x->edge);
      }
    };
    
    struct derivation_equal_type
    {
      bool operator()(const derivation_type* x, const derivation_type* y) const
      {
	return x->edge == y->edge && x->j == y->j;
      }
    };

#ifdef HAVE_TR1_UNORDERED_SET
    typedef std::tr1::unordered_set<const derivation_type*, derivation_hash_type, derivation_equal_type,
				    std::allocator<const derivation_type*> > derivation_set_unique_type;
#else
    typedef sgi::hash_set<const derivation_type*, derivation_hash_type, derivation_equal_type,
			  std::allocator<const derivation_type*> > derivation_set_unique_type;
#endif

    
    struct State
    {
      State() {}
      derivation_heap_type cand;
      derivation_list_type D;
      derivation_set_unique_type uniques;
    };
    
    typedef State state_type;
    typedef std::vector<state_type, std::allocator<state_type> > state_set_type;

    bool operator()(int v, int k, yield_type& yield)
    {
      const derivation_type* derivation = lazy_kth_best(v, k);
      if (derivation) {
	yield = derivation->yield;
	return true;
      } else {
	yield = yield();
	return false;
      }
    }
    
  private:
    
    const derivation_type* lazy_kth_best(int v, int k)
    {
      typedef std::vector<const yield_type*, std::allocator<const yield_type*> > yield_set_type;

      state_type & state = get_candidate(v);
      derivation_heap_type& cand = state.scand;
      derivation_list_type& D = state.D;

      yield_set_type yields;
      
      while (D.size() <= k) {
	
	if (D.size() > 0)
	  lazy_next(*D.back(), state);
	
	if (cand.size() > 0) {
	  std::pop_heap(cand.begin(), cand.end(), compare_heap_type());
	  const derivation_type* derivation = cand.back();
	  cand.pop_back();
	  
	  // perform traversal here...
	  
	  yields.clear();
	  for (int i = 0; i < derivation->edge->tail_ndoes.size(); ++ i)
	    yields.push_back(&lazy_kth_best(derivation->edge->tail_nodes[i], derivation->j[i])->yield);
	  
	  traversal(*(derivation->edge), &derivation->yield, yields.begin(), yields.end());
	  
	  // perform filtering here...!
	  // if we have duplicates, do not insert...
	  D.push_back(derivation);
	} else
	  break;
      }
      
      return (k < D.size() ? D[k] : 0);
    }
    
    void lazy_next(const derivation_type& derivation, state_type& state)
    {
      derivation_type query(derivation.j);
      index_set_type& j = query.j;
      
      for (int i = 0; i < j.size(); ++ i) {
	++ j[i];
	
	const derivation_type* antecedent = lazy_kth_best(derivation.edge->tail_nodes[i], j[i]);
	
	if (antecedent) {
	  query.edge = derivation.edge;
	  
	  if (state.uniques.find(&query) == state.uniques.end()) {
	    const derivation_type* derivation_new = make_derivation(*(derivation.edge), j);
	    
	    if (derivation_new) {
	      state.cand.push_back(derivation_new);
	      std::push_heap(state.cand.begin(), state.cand.end(), compare_heap_type());
	      state.uniques.insert(derivation_new);
	    }
	  }
	}
	-- j[i];
      }
    }

    const derivation_type* make_derivation(const edge_type& edge, const index_set_type& j)
    {
      derivations.push_back(derivation_type(j));
      
      derivation_type& derivation = derivations.back();
      
      derivation.score = function(edge);
      
      index_set_type::const_iterator iiter = j.begin();
      edge_type::node_set_type::const_iterator niter_end = edge.tail_nodes.end();
      for (edge_type::node_set_type::const_iterator niter = edge.tail_nodes.begin(); niter != niter_end; ++ niter, ++ iiter) {
	const derivation_type* antecedent = lazy_kth_best(*niter, *iiter);
	
	if (! antecedent) return 0;
	
	derivation.score *= antecedent->score;
      }
      
      return &derivation;
    }

    state_type& get_candidate(int v)
    {
      state_type& state = states[v];
      
      if (! state.D.empty() || ! state.cand.empty()) return state;
      
      const node_type& node = graph.nodes[v];
      
      node_type::edge_set_type::const_iterator eiter_end = node.in_edges.end();
      for (node_type::edge_set_type::const_iterator eiter = node.in_edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = graph.edges[*eiter];
	
	const index_set_type j(edge.tail_nodes.size(), 0);
	derivation_type* derivation = make_derivation(edge, j);
	
	if (! derivation)
	  throw std::runtime_error("no derivation?");
	
	state.cand.push_back(derivation);
      }
      
      // top k elements
      const size_type size = utils::bithack::min(k_prime, state.cand.size());
      
      std::nth_element(state.cand.begin(), state.cand.begin() + size, state.cand.end(), compare_derivation_type());
      state.cand.resize(size);
      
      // heapify
      std::make_heap(state.cand.begin(), state.cand.end(), compare_heap_type());
      
      return state;
    }
    
  private:
    const traversal_type traversal;
    const function_type  function;
    
    const hypergraph_type& graph;
    
    derivation_set_type derivations;
    state_set_type      states;
    
    const size_type k_prime;
  };
  
};

#endif
