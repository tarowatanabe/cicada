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
  // equality is defined by: state, source-edge-id, dot and j
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
    typedef hypergraph_type::rule_type rule_type;

    typedef hypergraph_type::feature_set_type feature_set_type;

    typedef Model model_type;
    
    typedef model_type::state_type     state_type;
    typedef model_type::state_set_type state_set_type;
    
    typedef Semiring semiring_type;
    typedef Semiring score_type;
    
    typedef Function function_type;
    
    struct Candidate
    {
      typedef Candidate candidate_type;
      
      state_type            state;
      const candidate_type* parent;
      
      const edge_type* in_edge;
      edge_type        out_edge;
      
      int dot;
      int dot_antecedent;
      
      score_type score;
      score_type estimate;
      
      Candidate() : parent(0), in_edge(0), dot(0), dot_antecedent(0) {}
      Candidate(const edge_type& __edge)
	: parent(0), in_edge(&__edge), out_edge(__edge), dot(0), dot_antecedent(0) {}
      Candidate(const candidate_type& __parent, const edge_type& __edge)
	: parent(&__parent), in_edge(&__edge), out_edge(__edge), dot(0), dot_antecedent(0) {}
    };
    
    typedef Candidate candidate_type;
    typedef utils::chunk_vector<candidate_type, 4096 / sizeof(candidate_type), std::allocator<candidate_type> > candidate_set_type;
    
    // hash and equal for keeping derivations
    struct candidate_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      candidate_hash_type(const model_type::state_hash& x) : hasher(x) {}
      
      model_type::state_hash hasher;

      size_t operator()(const candidate_type& x) const
      {
	return (x.in_edge == 0
		? size_t(0)
		: hasher_type::operator()(x.parent,
					  hasher_type::operator()(x.out_edge.tails.begin(), x.out_edge.tails.end(),
								  hasher_type::operator()(x.in_edge->id,
											  hasher_type::operator()(x.dot,
														  hasher(x.state))))));
      }
      
      size_t operator()(const candidate_type* x) const
      {
	return (x == 0 ? size_t(0) : operator()(*x));
      }
    };
    
    struct candidate_equal_type : public model_type::state_equal
    {
      candidate_equal_type(const model_type::state_equal& x) : model_type::state_equal(x) {}
      
      bool operator()(const candidate_type& x, const candidate_type& y) const
      {
	return x.in_edge == y.in_edge && (x.in_edge == 0
					  || (x.parent == y.parent
					      && x.dot == y.dot
					      && x.out_edge.tails == y.out_edge.tails
					      && model_type::state_equal::operator()(x.state, y.state)));
      }

      bool operator()(const candidate_type* x, const candidate_type* y) const
      {
	return (x == y) || (x && y && operator()(*x, *y));
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
    
    typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_base_type;
    typedef utils::b_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type, 512 / sizeof(const candidate_type*)> candidate_heap_type;
    
    typedef google::dense_hash_map<state_type, id_type, model_type::state_hash, model_type::state_equal > state_node_map_type;

    
    struct State
    {
      State(const size_type& hint, const size_type& state_size)
	: nodes(hint >> 1, model_type::state_hash(state_size), model_type::state_equal(state_size))
      {
	nodes.set_empty_key(state_type());
      }
      
      candidate_heap_type buf;
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
		
	states.clear();
	states.reserve(graph_in.nodes.size());
	states.resize(graph_in.nodes.size(), cand_state_type(pop_size_max >> 1, model.state_size()));
	
	// initialization...
	initialize_bins(graph_in);

	//std::cerr << "graph size: " << graph_in.nodes.size() << std::endl;
	
	for (int step = 0; step != static_cast<int>(graph_in.nodes.size()); ++ step)
	  process_bins(step, graph_in, graph_out);
	
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

	candidates.push_back(candidate_type(edge));

	candidate_type& candidate = candidates.back();
	
	feature_set_type estimates;
	model.apply_predict(candidate.state, node_states, candidate.out_edge, candidate.out_edge.features, estimates, true);
	
	candidate.score    = function(candidate.out_edge.features);
	candidate.estimate = function(estimates);
	candidate.estimate *= candidate.score;
	
	states.front().buf.push(&candidate);
      }
    }
    
    void process_bins(const int step, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      //std::cerr << "step: " << step << std::endl;

      for (size_type num_pop = 0; ! states[step].buf.empty() && num_pop != pop_size_max; ++ num_pop) {
	const candidate_type* item_top = states[step].buf.top();
	candidate_type* item = const_cast<candidate_type*>(item_top);
	
	states[step].buf.pop();
	
#if 0
	std::cerr << "head: " << item->in_edge->head
		  << " edge: " << item->in_edge->id
		  << " rule: " << *(item->in_edge->rule)
		  << " dot: " << item->dot
		  << std::endl;
#endif
	
	// we will iterate until completion...
	for (;;) {
	  const bool is_goal = (graph_in.goal == item->in_edge->head);
	  
 	  const rule_type::symbol_set_type& target = item->in_edge->rule->target;
	  
	  // scan... and score... state will be updated...
	  if (! target.empty() && item->dot < static_cast<int>(target.size()) && ! target[item->dot].is_non_terminal()) {
	    
	    if (item == item_top) {
	      candidates.push_back(*item_top);
	      item = &candidates.back();
	      item->state = model.clone(item->state);
	    }

	    feature_set_type features;
	    feature_set_type estimates;
	    
	    //std::cerr << "scan" << std::endl;
	    
	    model.apply_scan(item->state, node_states, item->out_edge, item->dot, features, estimates, is_goal);
	    
	    //std::cerr << "scan done" << std::endl;
	    
	    item->out_edge.features += features;
	    
	    const score_type score    = function(features);
	    const score_type estimate = function(estimates);
	    
	    item->score    *= score;
	    item->estimate *= score;
	    item->estimate *= estimate;
	    
	    // proceed dot(s)
	    for (/**/; item->dot < static_cast<int>(target.size()) && ! target[item->dot].is_non_terminal(); ++ item->dot);
	  }
	  
	  if (item->dot == static_cast<int>(target.size())) {
	    if (item == item_top) {
	      candidates.push_back(*item_top);
	      item = &candidates.back();
	      item->state = model.clone(item->state);
	    }
	    
	    // complete...
	    feature_set_type features;
	    feature_set_type estimates;
	    
	    // scoring
	    //std::cerr << "complete: " << (is_goal ? "goal" : "non-goal") << std::endl;
	    
	    model.apply_complete(item->state, node_states, item->out_edge, features, estimates, is_goal);

	    //std::cerr << "complete done" << std::endl;
	    
	    item->out_edge.features += features;
	    
	    const score_type score    = function(features);
	    const score_type estimate = function(estimates);
	    
	    item->score    *= score;
	    item->estimate *= score;
	    item->estimate *= estimate;
	    
	    //
	    // graph_out's node is differentiated by in_edge->head and current state...
	    //
	    
	    if (is_goal) {
	      if (! graph_out.is_valid()) {
		node_states.push_back(item->state);
		graph_out.goal = graph_out.add_node().id;
	      } else
		model.deallocate(item->state);
	      
	      edge_type& edge_new = graph_out.add_edge(item->out_edge);
	      graph_out.connect_edge(edge_new.id, graph_out.goal);
	      
	      // we will not use item any more...
	      // can we re-use this...?
	      
	      break;
	    } else {
	      
	      bool propagate = false;
	      
	      state_node_map_type::iterator siter = states[item->in_edge->head].nodes.find(item->state);
	      if (siter == states[item->in_edge->head].nodes.end()) {
		node_states.push_back(item->state);
		siter = states[item->in_edge->head].nodes.insert(std::make_pair(item->state, graph_out.add_node().id)).first;
		
		//std::cerr << "new node id: " << siter->second << std::endl;
		
		propagate = true;
	      } else 
		model.deallocate(item->state);
	      
	      edge_type& edge_new = graph_out.add_edge(item->out_edge);
	      graph_out.connect_edge(edge_new.id, siter->second);

	      if (! propagate) break;
	      
	      const score_type score    = item->score;
	      const score_type estimate = item->estimate;
	      
	      // some trick:
	      // item->state is either deleted or inserted in states[item->in_edge->head].nodes
	      // thus, we simply copy stat from item->parent...
	      // but reassigned from siter->first by cloning...
	      
	      *item = *(item->parent);
	      item->state = model.clone(siter->first);
	      item->score    *= score;
	      item->estimate *= estimate;
	      
#if 0
	      std::cerr << "parent head: " << item->in_edge->head
			<< " edge: " << item->in_edge->id
			<< " rule: " << *(item->in_edge->rule)
			<< " dot: " << item->dot
			<< std::endl;
#endif
	      
	      const rule_type::symbol_set_type& target = item->in_edge->rule->target;
	      
	      int antecedent_index = target[item->dot].non_terminal_index() - 1;
	      if (antecedent_index < 0)
		antecedent_index = item->dot_antecedent;
	      
	      item->out_edge.tails[antecedent_index] = siter->second;
	      
	      ++ item->dot;
	      ++ item->dot_antecedent;
	      
	      continue;
	    }
	    
	  } else {
	    // predict...
	    int antecedent_index = target[item->dot].non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      antecedent_index = item->dot_antecedent;
	    
	    const candidate_type& parent = *item;
	    const node_type& antecedent_node = graph_in.nodes[item->out_edge.tails[antecedent_index]];
	    
	    node_type::edge_set_type::const_iterator eiter_end = antecedent_node.edges.end();
	    for (node_type::edge_set_type::const_iterator eiter = antecedent_node.edges.begin(); eiter != eiter_end; ++ eiter) {
	      const edge_type& edge = graph_in.edges[*eiter];
	      
	      candidates.push_back(candidate_type(parent, edge));
	      
	      candidate_type& candidate = candidates.back();
	      
	      candidate.state = model.clone(item->state);
	      
	      feature_set_type estimates;
	      
	      model.apply_predict(candidate.state, node_states, candidate.out_edge, candidate.out_edge.features, estimates, false);
	      
	      const score_type score    = function(candidate.out_edge.features);
	      const score_type estimate = function(estimates);
	      
	      candidate.score    = parent.score;
	      candidate.estimate = parent.estimate;
	      
	      candidate.score    *= score;
	      candidate.estimate *= score;
	      candidate.estimate *= estimate;
	      
	      states[step + 1].buf.push(&candidate);
	    }
	    
	    break;
	  }
	}
      }
      
      // clear buf... at the same time, we will clear state associated with each candidate...
      
#if 0
      typename candidate_heap_type::const_iterator biter_end = states[step].buf.end();
      for (typename candidate_heap_type::const_iterator biter = states[step].buf.begin(); biter != biter_end; ++ biter) {
	model.deallocate((*biter)->state);
	
	// what to do with (*biter) == candidate_type* ??? 
	// I want to cache this for further reuse...
	// implement custom allocator similar to the state allocator in model?
	
      }
#endif
      
      // shrink wrap...
      states[step].buf.clear();
      candidate_heap_type(states[step].buf).swap(states[step].buf);
    };
    
  private:
    candidate_set_type  candidates;
    state_set_type      node_states;
    cand_state_set_type states;

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
