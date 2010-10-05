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
  //
  
  //
  // we need to keep state, source-edge-id, dot and j
  //
  
  // model.apply_predict   state-less feature application
  // model.apply_scan      state-full feature application with partial context
  // model.apply_complete  state-full feature application with full antecedent context
  
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
	return x.in_edge == y.in_edge && (x.in_edge == 0 || (x.dot == y.dot && x.out_edge.tails == y.out_edge.tails));
      }
    };
    typedef std::allocator<stack_state_type> stack_state_alloc_type;
    typedef utils::indexed_trie<stack_state_type, stack_state_hash_type, stack_state_equal_type, stack_state_alloc_type> stack_type;
    
    struct Candidate
    {
      state_type state;
      stack_type::id_type stack;
      
      feature_set_type features;
      
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
    };
    
    typedef State cand_state_type;
    typedef std::vector<cand_state_type, std::allocator<cand_state_type> > cand_state_set_type;
    
    ApplyIncremental(const model_type& _model,
		  const function_type& _function,
		  const int _pop_size_max)
      : model(_model),
	function(_function),
	pop_size_max(_pop_size_max)
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
	node_states.reserve(graph_in.nodes.size() * pop_size_max);

	node_states_coarse.clear();
	node_states_coarse.reserve(graph_in.nodes.size() * pop_size_max);

	node_maps.clear();
	node_maps.reserve(graph_in.nodes.size() * pop_size_max);
	
	states.clear();
	states.reserve(graph_in.nodes.size());
	states.resize(graph_in.nodes.size(), cand_state_type(pop_size_max >> 1, model.state_size()));
	
	// initialization...
	initialize_bins(graph_in);
	
	for (int cardinality = 0; cardinality != graph_in.nodes.size(); ++ cardinality)
	  process_bins(cardinality, graph_in, graph_out);
	
	//std::cerr << "topologically sort" << std::endl;
	
	// topologically sort...
	if (graph_out.is_valid())
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
	
	candidates.push_back(candidate_type(stack.push(stack_type::npos(), stack_state_type(edge))));
	
	candidate_type& candidate = candidates.back();
	
	// perform scoring by model.apply_predict()
	// new state is kept in candidate
	
	feature_set_type features(edge.features);
	feature_set_type estimates;
	model.apply_predict(candidate.state, node_states, stack[candidate.stack].out_edge, features, estimates, true);
	
	candidate.score    = function(features);
	candidate.estimate = function(estimates);
	candidate.estimate *= candidate.score;
	
	candidate.features.swap(features);
	
	states.front().buf.push(&candidate);
      }
    }
    
    void process_bins(const int cardinality, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      cand_state_type& cand_state = states[cardinality];
      
      for (size_type num_pop = 0; ! cand_state.buf.empty() && num_pop != pop_size_max; ++ num_pop) {
	candidate_type item = *(cand_state.buf.top());
	cand_state.buf.pop();

	const bool is_goal = (graph_in.goal == stack[item.stack].in_edge->head);
	
	// we will iterate until completion...
	for (;;) {
 	  const rule_type::symbol_set_type& target = stack[item.stack].in_edge->rule->target;
	  
	  int dot = stack[item.stack].dot;
	  int dot_antecedent = stack[item.stack].dot_antecedent;
	  
	  // scan... and score... state will be updated...
	  if (! target.empty() && target[dot].is_terminal()) {
	    feature_set_type features;
	    feature_set_type estimates;
	    
	    model.apply_scan(item.state, node_states, stack[item.stack].out_edge, dot, features, estimates, is_goal);
	    
	    const score_type score    = function(features);
	    const score_type estimate = function(estimates);
	    
	    item.features += features;
	    item.score    *= score;
	    item.estimate *= score;
	    item.estimate *= estimate;
	  }
	  
	  // proceed dot(s)
	  for (/**/; dot != target.size() && target[dot].is_terminal(); ++ dot);
	  
	  if (dot == target.size()) {
	    // complete...
	    
	    feature_set_type features;
	    feature_set_type estimates;
	    
	    // scoring
	    model.apply_complete(state, node_states, stack[item.stack].out_edge, features, estimates, is_goal);
	    
	    const score_type score    = function(features);
	    const score_type estimate = function(estimates);
	    
	    item.features += features;
	    item.score    *= score;
	    item.estimate *= score;
	    item.estimate *= estimate;
	    
	    //
	    // graph_out's node is differentiated by in_edge->head and current state...
	    //
	    
	    if (is_goal) {
	      if (! graph_out.is_valid()) {
		node_states.pupsh_back(state);
		graph_out.goal = graph_out.add_node().id;
	      } else
		model.deallocate(state);
	      
	      edge_tyep& edge_new = graph_out.add_edge(stack[item.stack].out_edge);
	      
	      edge_new.features.swap(item.features);
	      item.features.clear();
	      
	      graph_out.connect_edge(edge_new.id, graph_out.goal);
	      
	      break;
	    } else {
	      state_node_map_type::iterator siter = states[stack[item.stack].in_edge->head].nodes.find(state);
	      if (siter == states[stack[item.stack].in_edge->head].nodes.end()) {
		node_states.pupsh_back(state);
		siter = states[stack[item.stack].in_edge->head].nodes.insert(std::make_pair(state, graph_out.add_node().id)).first;
	      } else
		model.deallocate(state);
	      
	      edge_tyep& edge_new = graph_out.add_edge(stack[item.stack].out_edge);
	      
	      edge_new.features.swap(item.features);
	      item.features.clear();
	      
	      graph_out.connect_edge(edge_new.id, siter->second);
	      
	      // pop-stack...
	      const stack_type::id_type stack_parent = stack.parent(item.stack);
	      
	      stack_state_type stack_state(stack[stack_parent]);
	      
	      const rule_type::symbol_set_type& target = stack_state.in_edge->rule->target;
	      
	      int antecedent_index = target[stack_state.dot].non_terminal_index() - 1;
	      if (antecedent_index < 0)
		antecedent_index = stack_state.dot_antecedent;
	      
	      stack_state.out_edge.tails[antecedent_index] = siter->second;
	      
	      ++ stack_state.dot;
	      ++ stack_state.dot_ancedent;
	      
	      item.stack = stack.push(stack.parent(stack_parent), stack_state);
	      
	      continue;
	    }
	    
	  } else {
	    // predict...
	    int antecedent_index = target[dot].non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      antecedent_index = dot_antecedent;
	    
	    const node_type& antecedent_node = graph_in.nodes[stack[item.stack].out_edge.tails[antecedent_index]];
	    
	    node_type::edge_set_type::const_iterator eiter_end = antecedent_node.edges.end();
	    for (node_type::edge_set_type::const_iterator eiter = antecedent_node.edges.begin(); eiter != eiter_end; ++ eiter) {
	      const edge_type& edge = graph_in.edges[*eiter];
	      
	      candidates.push_back(candidate_type(stack.push(item.stack, stack_state_type(edge))));
	      candidate_type& candidate = candidates.back();
	      
	      feature_set_type features(edge.features);
	      feature_set_type estimates;
	      
	      candidate.state = model.clone(state);
	      
	      model.apply_predict(candidate.state, node_states, stack[candidate.stack].out_edges, features, estimates, false);
	      
	      const score_type score    = function(features);
	      const score_type estimate = function(estimates);
	      
	      candidate.score    = item.score;
	      candidate.estimate = item.estimate;
	      
	      candidate.score    *= score;
	      candidate.estimate *= score;
	      candidate.estimate *= estimate;
	      
	      candidate.features.swap(features);
	      
	      states[cardinality + 1].buf.push(&candidate);
	    }
	    
	    break;
	  }
	}
      }
      
      // clear buf... at the same time, we will clear state associated with each candidate...
      
      typename candidate_heap_type::const_iterator biter_end = state.buf.end();
      for (typename candidate_heap_type::const_iterator biter = state.buf.begin(); biter != biter_end; ++ biter) {
	model.deallocate((*biter)->state);
	
	// what to do with (*biter) == candidate_type* ??? 
	// I want to cache this for further reuse...
	// implement custom allocator similar to the state allocator in model?
	
      }
      
      // shrink wrap...
      state.buf.clear();
      candidate_heap_type(state.buf).swap(state.buf);
      
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
    size_type  pop_size_max;
  };
  
  template <typename Function>
  inline
  void apply_incremental(const Model& model, const HyperGraph& source, HyperGraph& target, const Function& func, const int pop_size)
  {
    ApplyIncremental<typename Function::value_type, Function>(model, func, pop_size)(source, target);
  }

  template <typename Function>
  inline
  void apply_incremental(const Model& model, HyperGraph& source, const Function& func, const int pop_size)
  {
    HyperGraph target;
    
    ApplyIncremental<typename Function::value_type, Function>(model, func, pop_size)(source, target);
    
    source.swap(target);
  }

};

#endif
