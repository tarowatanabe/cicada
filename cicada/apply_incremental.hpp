// -*- mode: c++ -*-

#ifndef __CICADA__APPLY_INCREMENTAL__HPP__
#define __CICADA__APPLY_INCREMENTAL__HPP__ 1

#include <cicada/apply_state_less.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/model.hpp>

#include <cicada/semiring/traits.hpp>

#include <google/dense_hash_set>

#include <utils/simple_vector.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/indexed_trie.hpp>

#include <utils/b_heap.hpp>

namespace cicada
{
  
  // implementation of incremental decoding, also knows as left-to-right
  // target generation
  
  //
  // we need feature application in predict/scan/complete
  //
  // we will use beam to hold all the expanded candidates, which is organized in vector (or use unique merging...?)
  // 
  // each candidate hold partially instantiated edge (in target-side left-to-right ordering)
  //   note: do we allow source-yield incrementatl parsing?
  //
  // How to represent "stack"?
  // indexed-trie of pointer of candidate
  // equality is defined by: source-edge-id, dot, j
  //
  
  template <typename Semiring, typename Function>
  struct ApplyIncremental
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

    struct StackState
    {
      const edge_type* in_edge;
      edge_type        out_edge;
      
      int dot;
      int dot_antecedent;

      StackState() : in_edge(0), dot(0), dot_antecedent(0) {}
      StackState(const edge_type& __edge)
	: in_edge(&__edge), out_edge(__edge), dot(0), dot_antecedent(0) {}
    };
    
    typedef StackState stack_state_type;
    
    struct stack_state_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const stack_state_type& x) const
      {
	return (x.in_edge == 0
		? size_t(0)
		: hasher_type::operator()(x.out_edge.tails.begin(), x.out_edge.tails.end(),
					  hasher_type::operator()(x.in_edge->id, x.dot)));
      }
    };
    struct stack_state_equal_type
    {
      bool operator()(const stack_state_type& x, const stack_state_type& y) const
      {
	return x.in_edge == y.in_edge && x.in_edge && y.in_edge && x.dot == y.dot && x.out_edge.tails == y.out_edge.tails;
      }
    };
    typedef std::allocator<stack_state_type> stack_state_alloc_type;
    typedef utils::indexed_trie<stack_state_type, stack_state_hash_type, stack_state_equal_type, stack_state_alloc_type> stack_type;
    
    struct Candidate
    {
      id_type node;
      
      state_type state;
      stack_type::id_type stack;
      
      score_type score;
      score_type estimate;
      
      Candidate() : stack(stack_type::npos()) {}
      Candidate(const stack_type::id_type& __stack) : stack(__stack) {}
    };
    
    typedef Candidate candidate_type;
    typedef utils::chunk_vector<candidate_type, 4096 / sizeof(candidate_type), std::allocator<candidate_type> > candidate_set_type;
    
    // hash and equal for keeping derivations
    struct candidate_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      model_type::state_hash state_hash;
      
      candidate_hash_type(const model_type::state_hash& __state_hash) : state_hash(__state_hash) {}
      
      size_t operator()(const candidate_type* x) const
      {
	return (x == 0 ? size_t(0) : hasher_type::opeator()(x->stack, state_hash(x->state)));
      }
    };
    
    struct candidate_equal_type : public model_type::state_equal
    {
      candidate_equal_type(const model_type::state_equal& x) : model_type::state_equal(x) {}
      
      bool operator()(const candidate_type* x, const candidate_type* y) const
      {
	return (x == y) || (x && y && x->stack == y->stack && model_type::state_equal::operator()(x->state, y->state));
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
    
    typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_list_type;
    
    typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_base_type;
    typedef utils::b_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type, 512 / sizeof(const candidate_type*)> candidate_heap_type;
    
    typedef google::dense_hash_map<state_type, id_type, model_type::state_hash, model_type::state_equal > state_node_map_type;
    typedef google::dense_hash_set<const candidate_type*, candidate_hash_type, candidate_equal_type > candidate_set_unique_type;

    typedef std::vector<id_type, std::allocator<id_type> > node_map_type;

    struct State
    {
      State(const size_type& hint, const size_type& state_size)
	: nodes(hint >> 1, model_type::state_hash(state_size), model_type::state_equal(state_size)),
	  nodes_coarse(hint, model_type::state_hash(state_size), model_type::state_equal(state_size)),
	  fired(false)
      {
	nodes.set_empty_key(state_type());
	nodes_coarse.set_empty_key(state_type());
	
	uniques.set_empty_key(0);
      }
      
      candidate_heap_type cand;
      candidate_heap_type buf;
      
      candidate_list_type D;
      candidate_set_unique_type uniques;
      
      state_node_map_type nodes;
      state_node_map_type nodes_coarse;

      bool fired;
    };
    
    typedef State cand_state_type;
    typedef std::vector<cand_state_type, std::allocator<cand_state_type> > cand_state_set_type;
    
    ApplyIncremental(const model_type& _model,
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

	node_maps.clear();
	node_maps.reserve(graph_in.nodes.size() * cube_size_max);
	
	states.clear();
	states.reserve(graph_in.nodes.size());
	states.resize(graph_in.nodes.size(), cand_state_type(cube_size_max >> 1, model.state_size()));
	
	// initialization...
	initialize_bins(graph_in);
	
	for (int cardinality = 0; cardinality != graph_in.nodes.size(); ++ cardinality)
	  process_bins(cardinality, graph_in, graph_out);
	
	//std::cerr << "topologically sort" << std::endl;
	
	// topologically sort...
	graph_out.topologically_sort();
	
	//std::cerr << "topologically sort: end" << std::endl;
	
	// re-initialize again...
	const_cast<model_type&>(model).initialize();
      }
    };
    
  private:

    void initialize_bins(const hypergraph_type& graph)
    {
      node_type::edge_set_type::const_iterator eiter_end = graph.nodes[graph.goal].edges.end();
      for (node_type::edge_set_type::const_iterator eiter = graph.nodes[graph.goal].edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = graph.edges[*eiter];
	
	const stack_type::id_type stack_id = stack.push(stack_type::npos(), stack_state_type(*eiter, 0));
	
	candidates.push_back(candidate_type(edge, stack_id, index_set_type(edge.tails.begin(), edge.tails.end())));
	
	candidate_type& candidate = candidates.back();
	
	// perform scoring by model.apply_predict()
	// state is kept in candidate
	
	candidate.stack = stack_id;
	
	states.front().buf.push(&candidate);
      }
    }
    
    void process_bins(const int cardinality, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      cand_state_type& state = states[cardinality];
      
      for (size_type num_pop = 0; ! state.buf.empty() && num_pop != pop_size_max; ++ num_pop) {
	const candidate_type* item = state.buf.top();
	state.buf.pop();
	
	// we will iterate until completion...
	for (;;) {
	  const rule_type::symbol_set_type& target = item->in_edge->rule->target;
	  
	  int dot = item->dot;
	  int dot_antecedent = stack[item->stack].second;
	  
	  // scan... and score...
	  if (target[dot].is_terminal())
	    model.apply_scan();
	  
	  for (/**/; dot != target.size() && target[dot].is_terminal(); ++ dot);
	  
	  if (dot == target.size()) {
	    // complete...
	    // pop-stack and loop again...
	    
	    continue;
	  } else {
	    // predict...
	    
	    break;
	  }
	}
	
      }
      
      // clear buf... at the same time, we will clear state associated with each candidate...
      
      typename candidate_heap_type::const_iterator biter_end = state.buf.end();
      for (typename candidate_heap_type::const_iterator biter = state.buf.begin(); biter != biter_end; ++ biter)
	model.deallocate((*biter)->state);
      
      // shrink wrap...
      state.buf.clear();
      candidate_heap_type(state.buf).swap(state.buf);
      
      
      for () {
	
	const candidate_type* item = *biter;
	
	// we will iterate completion...
	for (;;) {
	  const rule_type::symbol_set_type& target = item->in_edge->rule->target;
	  
	  int dot = item->dot;
	  
	  feature_set_type estimates;
	  
	  score_type score    = semiring::traits<score_type>::one();
	  score_type estimate = semiring::traits<score_type>::one();
	  
	  // scan
	  if (target[dot + 1].is_terminal())
	    model.apply_scan();
	  
	  for (/**/; dot != target.size() && target[dot + 1].is_terminal(); ++ dot);
	  
	  
	  if (dot == target.size()) {
	    // complete
	    
	    // pop-stack and loop again...
	    
	    continue;
	  } else {
	    // predict
	    // TODO: pop from stack and proceed "dot"
	    item = item->stack;
	    ++ item->dot;
	  
	    break;
	  }
	}
      }
      
      
    };
    
    void lazy_jth_best(id_type v, int j, const hypergraph_type& graph, hypergraph_type& graph_out)
    {
      //std::cerr << "node: " << v << std::endl;
      
      const bool is_goal = graph.goal == v;
      
      cand_state_type& state = states[v];
      
      if (! state.fired) {
	const node_type& node = graph.nodes[v];
	
	// for each edge in v
	//   fire(edge, 0, cand)
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const edge_type& edge = graph.edges[*eiter];
	  const index_set_type j(edge.tails.size(), 0);
	  
	  fire(edge, j, state, graph, graph_out);
	}
	
	state.fired = true;
      }
      
      while (state.D.size() <= j && state.buf.size() + state.D.size() < cube_size_max && ! state.cand.empty()) {
	// pop-min
	const candidate_type* item = state.cand.top();
	state.cand.pop();
	
	// push item to buffer
	push_buf(*item, state, is_goal, graph_out);
	
	// push succ
	push_succ(*item, state, graph, graph_out);
	
	// enum item with current bound
	if (! state.cand.empty())
	  enum_item(state, state.cand.top()->estimate, is_goal, graph_out);
      }
      
      // enum item with zero bound
      enum_item(state, semiring::traits<score_type>::zero(), is_goal, graph_out);
    }

    void fire(const edge_type& edge, const index_set_type& j, cand_state_type& state, const hypergraph_type& graph, hypergraph_type& graph_out)
    {
      candidate_type query(j);
      query.in_edge = &edge;
      
      if (state.uniques.find(&query) != state.uniques.end()) return;
      
      // for each edge, 
      index_set_type::const_iterator iiter = j.begin();
      edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
      for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter, ++ iiter) {
	lazy_jth_best(*niter, *iiter, graph, graph_out);
	
	if (states[*niter].D.size() <= *iiter) return;
      }
      
      // push cand
      const candidate_type* item = make_candidate(edge, j, state, graph.goal == edge.head, graph_out);
      
      state.uniques.insert(item);
      state.cand.push(item);
    }
    
    void push_succ(const candidate_type& item, cand_state_type& state, const hypergraph_type& graph, hypergraph_type& graph_out)
    {
      index_set_type j = item.j;
      
      // for each i in 1 ... |e|
      //   fire(e, j + b_i, cand)
      for (int i = 0; i < item.j.size(); ++ i) {
	const int j_i_prev = j[i];
	++ j[i];
	
	fire(*item.in_edge, j, state, graph, graph_out);
	
	j[i] = j_i_prev;
      }
    }
    
    void enum_item(cand_state_type& state, const score_type bound, const bool is_goal, hypergraph_type& graph_out)
    {
      // while |buf| and min(buf) < bound (min-cost)
      //  append pop-min to D
      while (! state.buf.empty() && state.buf.top()->estimate > bound) {
	const candidate_type* item = state.buf.top();
	state.buf.pop();
	
	candidate_type& candidate = const_cast<candidate_type&>(*item);
	
	// we will add new node/edge here 
	// If possible, state merge
	if (is_goal) {
	  if (! graph_out.is_valid()) {
	    node_maps.push_back(candidate.node);
	    node_states.push_back(candidate.state);
	    
	    graph_out.goal = graph_out.add_node().id;
	  } else
	    model.deallocate(candidate.state);
	  
	  candidate.node = graph_out.goal;
	} else {
	  // we will merge states, but do not merge score/estimates, since we
	  // are enumerating jthe best derivations... is this correct?
	  state_node_map_type::iterator siter = state.nodes.find(candidate.state);
	  if (siter == state.nodes.end()) {
	    node_maps.push_back(candidate.node);
	    node_states.push_back(candidate.state);
	    
	    siter = state.nodes.insert(std::make_pair(candidate.state, graph_out.add_node().id)).first;
	  } else
	    model.deallocate(candidate.state);
	  
	  candidate.node = siter->second;
	}
	
	edge_type& edge = graph_out.add_edge(candidate.out_edge);
	graph_out.connect_edge(edge.id, candidate.node);
	
	state.D.push_back(item);
      }
    }

    void push_buf(const candidate_type& __item, cand_state_type& state, const bool is_goal, hypergraph_type& graph_out)
    {
      // push this item into state.buf with "correct" score

      candidate_type& candidate = const_cast<candidate_type&>(__item);
      
      candidate.score = semiring::traits<score_type>::one();
      candidate.estimate = semiring::traits<score_type>::one();
      for (int i = 0; i < candidate.j.size(); ++ i) {
	const candidate_type& antecedent = *states[candidate.in_edge->tails[i]].D[candidate.j[i]];
	
	// assign real-node-id!
	candidate.out_edge.tails[i] = antecedent.node;
	candidate.score *= antecedent.score;
      }
      
      feature_set_type estimates;
      candidate.state = model.apply(node_states, candidate.out_edge, estimates, is_goal);
      
      candidate.score    *= function(candidate.out_edge.features);
      candidate.estimate *= function(estimates);
      candidate.estimate *= candidate.score;
      
      state.buf.push(&candidate);
    }
    
    const candidate_type* make_candidate(const edge_type& edge, const index_set_type& j, cand_state_type& state, const bool is_goal, hypergraph_type& graph_out)
    {
      candidates.push_back(candidate_type(edge, j));
      
      candidate_type& candidate = candidates.back();
      
      candidate.out_edge.tails = edge_type::node_set_type(j.size());
      
      candidate.score = semiring::traits<score_type>::one();
      candidate.estimate = semiring::traits<score_type>::one();
      for (int i = 0; i < j.size(); ++ i) {
	const candidate_type& antecedent = *states[edge.tails[i]].D[j[i]];
	
	// assign coarse node id
	candidate.out_edge.tails[i] = node_maps[antecedent.node];
	candidate.score *= antecedent.score;
      }
      
      // perform "estimated" coarse model application
      feature_set_type estimates;
      const state_type node_state = model.apply_coarse(node_states_coarse, candidate.out_edge, estimates, is_goal);
      
      candidate.score    *= function(candidate.out_edge.features);
      candidate.estimate *= function(estimates);
      candidate.estimate *= candidate.score;
      
      // no state merging...
      //candidate.node = node_states_coarse.size();
      //node_states_coarse.push_back(node_state);
      
      // state merging... so that we may reuse state structure
      state_node_map_type::iterator siter = state.nodes_coarse.find(node_state);
      if (siter == state.nodes_coarse.end()) {
	siter = state.nodes_coarse.insert(std::make_pair(node_state, node_states_coarse.size())).first;
	node_states_coarse.push_back(node_state);
      } else
	model.deallocate(node_state);
      
      candidate.node = siter->second;
      
      return &candidate;
    };
    
    
  private:
    candidate_set_type  candidates;
    state_set_type      node_states;
    state_set_type      node_states_coarse;
    cand_state_set_type states;

    node_map_type       node_maps;

    stack_type stack;

    const model_type& model;
    const function_type& function;
    size_type  cube_size_max;
  };
  
  template <typename Function>
  inline
  void apply_incremental(const Model& model, const HyperGraph& source, HyperGraph& target, const Function& func, const int cube_size)
  {
    ApplyIncremental<typename Function::value_type, Function>(model, func, cube_size)(source, target);
  }

  template <typename Function>
  inline
  void apply_incremental(const Model& model, HyperGraph& source, const Function& func, const int cube_size)
  {
    HyperGraph target;
    
    ApplyIncremental<typename Function::value_type, Function>(model, func, cube_size)(source, target);
    
    source.swap(target);
  }

};

#endif
