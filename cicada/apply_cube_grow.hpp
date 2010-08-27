// -*- mode: c++ -*-

#ifndef __CICADA__APPLY_CUBE_GROW__HPP__
#define __CICADA__APPLY_CUBE_GROW__HPP__ 1

#include <cicada/apply_state_less.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/model.hpp>

#include <cicada/semiring/traits.hpp>

#include <google/dense_hash_set>

#include <utils/simple_vector.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>

#include <utils/b_heap.hpp>

namespace cicada
{
  
  // faith full implementation of 
  //
  // @InProceedings{huang-chiang:2007:ACLMain,
  //  author    = {Huang, Liang  and  Chiang, David},
 //  title     = {Forest Rescoring: Faster Decoding with Integrated Language Models},
  //  booktitle = {Proceedings of the 45th Annual Meeting of the Association of Computational Linguistics},
  //  month     = {June},
  //  year      = {2007},
  //  address   = {Prague, Czech Republic},
  //  publisher = {Association for Computational Linguistics},
  //  pages     = {144--151},
  //  url       = {http://www.aclweb.org/anthology/P07-1019}
  //  }
  //

  
  // semiring and function to compute semiring from a feature vector

  template <typename Semiring, typename Function>
  struct ApplyCubeGrow
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type   id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    typedef hypergraph_type::feature_set_type feature_set_type;

    typedef Model model_type;
    
    typedef model_type::state_type     state_type;
    typedef model_type::state_set_type state_set_type;
    
    typedef Semiring semiring_type;
    typedef Semiring score_type;
    
    typedef Function function_type;
    
    typedef utils::simple_vector<int, std::allocator<int> > index_set_type;
    
    struct Candidate
    {
      id_type node;
      
      const edge_type* in_edge;
      edge_type        out_edge;
      
      index_set_type j;
      
      score_type score;
      score_type estimate;
      
      Candidate(const index_set_type& __j)
	: in_edge(0), j(__j) {}

      Candidate(const edge_type& __edge, const index_set_type& __j)
	: in_edge(&__edge), out_edge(__edge), j(__j) {}
    };

    typedef Candidate candidate_type;
    typedef utils::chunk_vector<candidate_type, 4096 / sizeof(candidate_type), std::allocator<candidate_type> > candidate_set_type;
        
    
    struct candidate_hash_type : public utils::hashmurmur<size_t>
    {
      size_t operator()(const candidate_type* x) const
      {
	return (x == 0 ? size_t(0) : utils::hashmurmur<size_t>::operator()(x->j.begin(), x->j.end(), x->in_edge->id));
      }
    };
    struct candidate_equal_type
    {
      bool operator()(const candidate_type* x, const candidate_type* y) const
      {
	return (x == y) || (x && y && x->in_edge->id == y->in_edge->id && x->j == y->j);
      }
    };
    
    struct compare_heap_type
    {
      // we use less, so that when popped from heap, we will grab "greater" in back...
      bool operator()(const candidate_type* x, const candidate_type* y) const
      {
	return x->estimate < y->estimate;
      }
    };

    struct node_score_type
    {
      id_type node;
      score_type score;
      score_type estimate;

      node_score_type() : node(), score(), estimate() {}
      
      node_score_type(const id_type __node, const score_type& __score, const score_type& __estimate)
	: node(__node), score(__score), estimate(__estimate) {}
    };
    
    typedef std::vector<node_score_type, std::allocator<node_score_type> > node_score_list_type;
    
    typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_base_type;
    typedef utils::b_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type, 512 / sizeof(const candidate_type*)> candidate_heap_type;
    
    typedef google::dense_hash_set<const candidate_type*, candidate_hash_type, candidate_equal_type > candidate_set_unique_type;

    struct State
    {
      State() : fired(false) { uniques.set_empty_key(0); }
      
      candidate_heap_type cand;
      candidate_heap_type buf;
      
      node_score_list_type D;
      candidate_set_unique_type uniques;

      bool fired;
    };
    
    typedef State cand_state_type;
    typedef std::vector<cand_state_type, std::allocator<cand_state_type> > cand_state_set_type;
    
    ApplyCubeGrow(const model_type& _model,
		  const function_type& _function,
		  const int _cube_size_max)
      : model(_model),
	function(_function),
	cube_size_max(_cube_size_max)
    { 
      
    }
    
    void operator()(const hypergraph_type& graph_in,
		    hypergraph_type&       graph_out)
    {
      const_cast<model_type&>(model).initialize();

      graph_out.clear();

      if (model.is_stateless()) {
	ApplyStateLess __applier(model);
	__applier(graph_in, graph_out);
      } else {
	candidates.clear();
	
	node_states.clear();
	node_states.reserve(graph_in.nodes.size() * cube_size_max);

	node_states_coarse.clear();
	node_states_coarse.reserve(graph_in.nodes.size() * cube_size_max);
	
	states.clear();
	states.reserve(graph_in.nodes.size());
	states.resize(graph_in.nodes.size());
	
	for (int j = 0; j < cube_size_max; ++ j) {
	  const size_type edge_size = graph_out.edges.size();
	  lazy_jth_best(graph_in.goal, j, graph_in, graph_out);

	  // no new edges inserted...
	  if (graph_out.edges.size() == edge_size) break;
	}

	//std::cerr << "topologically sort" << std::endl;
	
	// topologically sort...
	graph_out.topologically_sort();
	
	//std::cerr << "topologically sort: end" << std::endl;
	
	// re-initialize again...
	const_cast<model_type&>(model).initialize();
      }
    };
    
  private:
    
    void lazy_jth_best(id_type v, int j, const hypergraph_type& graph, hypergraph_type& graph_out)
    {
      //std::cerr << "node: " << v << std::endl;
      
      cand_state_type& state = states[v];
      
      if (! state.fired) {
	//std::cerr << "initial fire: " << v << std::endl;
	
	const node_type& node = graph.nodes[v];
	
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const edge_type& edge = graph.edges[*eiter];
	  const index_set_type j(edge.tails.size(), 0);
	  
	  fire(edge, j, state, graph, graph_out);
	}
	
	state.fired = true;
      }
      
      //std::cerr << "enumerate: " << v << std::endl;
      
      while (state.D.size() <= j && state.buf.size() + state.D.size() < cube_size_max && ! state.cand.empty()) {
	const candidate_type* item = state.cand.top();
	state.cand.pop();
	
	push_buf(*item, state, graph.goal == v, graph_out);
	
	push_succ(*item, state, graph, graph_out);
	
	if (! state.cand.empty())
	  enum_item(state, state.cand.top()->estimate, graph_out);
      }
      
      enum_item(state, semiring::traits<score_type>::zero(), graph_out);
    }

    void fire(const edge_type& edge, const index_set_type& j, cand_state_type& state, const hypergraph_type& graph, hypergraph_type& graph_out)
    {
      //std::cerr << "fire for: " << (*edge.rule) << std::endl;
      
      candidate_type query(j);
      query.in_edge = &edge;
      
      if (state.uniques.find(&query) != state.uniques.end()) return;
      
      index_set_type::const_iterator iiter = j.begin();
      edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
      for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter, ++ iiter) {
	lazy_jth_best(*niter, *iiter, graph, graph_out);
	
	if (states[*niter].D.size() <= *iiter) return;
      }
      
      //std::cerr << "fired: " << (*edge.rule) << std::endl;
      
      const candidate_type* item = make_candidate(edge, j, graph.goal == edge.head, graph_out);
      
      state.uniques.insert(item);
      state.cand.push(item);
    }
    
    void push_succ(const candidate_type& item, cand_state_type& state, const hypergraph_type& graph, hypergraph_type& graph_out)
    {
      //std::cerr << "push succ: " << (*item.in_edge->rule) << std::endl;
      
      index_set_type j = item.j;
      
      for (int i = 0; i < item.j.size(); ++ i) {
	const int j_i_prev = j[i];
	++ j[i];
	
	fire(*item.in_edge, j, state, graph, graph_out);
	
	j[i] = j_i_prev;
      }
    }
    
    void enum_item(cand_state_type& state, const score_type bound, hypergraph_type& graph_out)
    {
      while (! state.buf.empty() && state.buf.top()->estimate > bound) {
	const candidate_type* item = state.buf.top();
	state.buf.pop();
	
	state.D.push_back(node_score_type(item->node, item->score, item->estimate));
	
	edge_type& edge = graph_out.add_edge(item->out_edge);
	graph_out.connect_edge(edge.id, item->node);
      }
    }

    void push_buf(const candidate_type& __item, cand_state_type& state, const bool is_goal, hypergraph_type& graph_out)
    {
      // push this item into state.buf with "correct" score

      //std::cerr << "push buf for: " << *(__item.in_edge->rule) << std::endl;

      candidate_type& candidate = const_cast<candidate_type&>(__item);
      
      candidate.score = semiring::traits<score_type>::one();
      candidate.estimate = semiring::traits<score_type>::one();
      for (int i = 0; i < candidate.j.size(); ++ i) {
	const node_score_type& antecedent = states[candidate.in_edge->tails[i]].D[candidate.j[i]];
	
	candidate.score *= antecedent.score;
      }

      //std::cerr << "apply: " << std::endl;
      
      feature_set_type estimates;
      const state_type node_state = model.apply(node_states, candidate.out_edge, estimates, is_goal);
      
      candidate.score    *= function(candidate.out_edge.features);
      candidate.estimate *= function(estimates);
      candidate.estimate *= candidate.score;

      if (candidate.node >= node_states.size())
	node_states.resize(candidate.node + 1);
      
      if (node_states[candidate.node].empty())
	node_states[candidate.node] = node_state;
      else
	model.deallocate(node_state);
      
      state.buf.push(&candidate);

      //std::cerr << "end push buf" << std::endl;
    }
    
    const candidate_type* make_candidate(const edge_type& edge, const index_set_type& j, const bool is_goal, hypergraph_type& graph_out)
    {
      //std::cerr << "make candidate for: " << *(edge.rule) << std::endl;
      
      candidates.push_back(candidate_type(edge, j));
      
      candidate_type& candidate = candidates.back();
      
      candidate.out_edge.tails = edge_type::node_set_type(j.size());
      
      candidate.score = semiring::traits<score_type>::one();
      candidate.estimate = semiring::traits<score_type>::one();
      for (int i = 0; i < j.size(); ++ i) {
	const node_score_type& antecedent = states[edge.tails[i]].D[j[i]];
	
	candidate.out_edge.tails[i] = antecedent.node;
	candidate.score *= antecedent.score;
      }
      
      // perform "estimated" model application

      feature_set_type estimates;
      const state_type node_state = model.apply_coarse(node_states_coarse, candidate.out_edge, estimates, is_goal);

      candidate.score    *= function(candidate.out_edge.features);
      candidate.estimate *= function(estimates);
      candidate.estimate *= candidate.score;
      
      if (is_goal) {
	if (! graph_out.is_valid()) {
	  graph_out.goal = graph_out.add_node().id;
	  node_states_coarse.push_back(node_state);
	} else
	  model.deallocate(node_state);
	
	candidate.node = graph_out.goal;
      } else {
	candidate.node = graph_out.add_node().id;
	node_states_coarse.push_back(node_state);
      }
      
      return &candidate;
    };
    
    
  private:
    candidate_set_type  candidates;
    state_set_type      node_states;
    state_set_type      node_states_coarse;
    cand_state_set_type states;
    
    const model_type& model;
    const function_type& function;
    size_type  cube_size_max;
  };
  
  template <typename Function>
  inline
  void apply_cube_grow(const Model& model, const HyperGraph& source, HyperGraph& target, const Function& func, const int cube_size)
  {
    ApplyCubeGrow<typename Function::value_type, Function>(model, func, cube_size)(source, target);
  }

  template <typename Function>
  inline
  void apply_cube_grow(const Model& model, HyperGraph& source, const Function& func, const int cube_size)
  {
    HyperGraph target;
    
    ApplyCubeGrow<typename Function::value_type, Function>(model, func, cube_size)(source, target);
    
    source.swap(target);
  }

};

#endif
