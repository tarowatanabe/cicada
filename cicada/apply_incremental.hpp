// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__APPLY_INCREMENTAL__HPP__
#define __CICADA__APPLY_INCREMENTAL__HPP__ 1

#include <cicada/apply_state_less.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/model.hpp>
#include <cicada/inside_outside.hpp>

#include <cicada/semiring/traits.hpp>
#include <cicada/semiring/tropical.hpp>

#include <utils/dense_hash_map.hpp>
#include <utils/dense_hash_set.hpp>
#include <utils/simple_vector.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/indexed_trie.hpp>
#include <utils/b_heap.hpp>
#include <utils/std_heap.hpp>
#include <utils/bithack.hpp>

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

    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;
    
    typedef feature_set_type::feature_type     feature_type;
    typedef attribute_set_type::attribute_type attribute_type;
    
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
      
      Candidate()
	: state(), parent(0), in_edge(0), dot(0), dot_antecedent(0) {}
      Candidate(const edge_type& __edge)
	: state(), parent(0), in_edge(&__edge), out_edge(__edge), dot(0), dot_antecedent(0) {}
      Candidate(const candidate_type& __parent, const edge_type& __edge)
	: state(), parent(&__parent), in_edge(&__edge), out_edge(__edge), dot(0), dot_antecedent(0) {}
    };
    
    typedef Candidate candidate_type;
    typedef utils::chunk_vector<candidate_type, 4096 / sizeof(candidate_type), std::allocator<candidate_type> > candidate_set_type;
    
    // hash and equal for keeping derivations
    // Note: the state is omitted... but this is implicitly encoded by tails and dot position...
    struct candidate_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      size_t operator()(const candidate_type& x) const
      {
	return hasher_type::operator()(x.parent,
				       hasher_type::operator()(x.out_edge.tails.begin(), x.out_edge.tails.end(),
							       hasher_type::operator()(x.in_edge, x.dot)));
      }
      
      size_t operator()(const candidate_type* x) const
      {
	return (x == 0 ? size_t(0) : operator()(*x));
      }
    };
    
    struct candidate_equal_type
    {
      bool operator()(const candidate_type& x, const candidate_type& y) const
      {
	return (x.in_edge == y.in_edge
		&& x.parent == y.parent
		&& x.dot == y.dot
		&& x.out_edge.tails == y.out_edge.tails);
      }

      bool operator()(const candidate_type* x, const candidate_type* y) const
      {
	return (x == y) || (x && y && operator()(*x, *y));
      }
    };
    
    typedef typename utils::dense_hash_set<const candidate_type*, candidate_hash_type, candidate_equal_type>::type candidate_unique_set_type;
    
    struct compare_heap_type
    {
      // we use less, so that when popped from heap, we will grab "greater" in back...
      bool operator()(const candidate_type* x, const candidate_type* y) const
      {
	return x->score < y->score;
      }
    };
    
    typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_base_type;
    //typedef utils::b_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type, 512 / sizeof(const candidate_type*)> candidate_heap_type;
    typedef utils::std_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type> candidate_heap_type;
    typedef utils::chunk_vector<candidate_heap_type, 4096/sizeof(candidate_heap_type), std::allocator<candidate_heap_type> > candidate_heap_set_type;
    
    typedef std::pair<const candidate_type*, state_type> stack_state_type;
    
    struct stack_state_hash_type : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;

      model_type::state_hash hasher;
      
      stack_state_hash_type(const size_type& state_size)
	: hasher(state_size) {}
      
      size_type operator()(const stack_state_type& x) const
      {
	return hasher_type::operator()(x.first, hasher(x.second));
      }
    };
    
    struct stack_state_equal_type : public model_type::state_equal
    {
      stack_state_equal_type(const size_type& state_size)
	: model_type::state_equal(state_size) {}
      
      bool operator()(const stack_state_type& x, const stack_state_type& y) const
      {
	return x.first == y.first && model_type::state_equal::operator()(x.second, y.second);
      }
    };
    
    typedef typename utils::dense_hash_map<stack_state_type, id_type, stack_state_hash_type, stack_state_equal_type>::type state_node_map_type;
    
    struct State : public state_node_map_type
    {
      State(const size_type& hint, const size_type& state_size)
	: state_node_map_type(hint >> 1, stack_state_hash_type(state_size), stack_state_equal_type(state_size))
      {
	state_node_map_type::set_empty_key(stack_state_type(0, state_type()));
      }
    };
    
    typedef State cand_state_type;
    typedef std::vector<cand_state_type, std::allocator<cand_state_type> > cand_state_set_type;    
    
    struct max_edge_function
    {
      typedef cicada::semiring::Tropical<int> value_type;
      
      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	return cicada::semiring::traits<value_type>::exp(1);
      }
    };
    typedef typename max_edge_function::value_type count_type;
    
    typedef std::vector<count_type, std::allocator<count_type> > count_set_type;
    
    struct __attribute_integer : public boost::static_visitor<cicada::AttributeVector::int_type>
    {
      attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
      attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return 0; }
      attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return 0; }
    };

    ApplyIncremental(const model_type& _model,
		     const function_type& _function,
		     const int _pop_size_max)
      : model(_model),
	function(_function),
	pop_size_max(_pop_size_max)
	//attr_scan("incremental-scan"),
	//attr_complete("incremental-complete"),
	//attr_predict("incremental-predict")
    { 
      predictions.set_empty_key(0);
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
	candidate_list = 0;
	
	node_states.clear();
	node_states.reserve(graph_in.nodes.size() * pop_size_max);

	states.clear();
	states.reserve(graph_in.nodes.size());
	states.resize(graph_in.nodes.size(), cand_state_type(pop_size_max >> 1, model.state_size()));
	
	buf.clear();
	buf.resize(graph_in.nodes.size());
	
	counts_inside.clear();
	counts_inside.reserve(graph_in.nodes.size());
	counts_inside.resize(graph_in.nodes.size());

	//counts_outside.clear();
	//counts_outside.reserve(graph_in.nodes.size());
	//counts_outside.resize(graph_in.nodes.size());
	
	cicada::inside(graph_in, counts_inside, max_edge_function());
	//cicada::outside(graph_in, counts_inside, counts_outside, max_edge_function());
	
	// initialization...
	initialize_bins(graph_in);
	
	//std::cerr << "graph size: " << graph_in.nodes.size() << std::endl;
	
	for (size_t step = 0; step != buf.size(); ++ step)
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

	candidate_type& candidate = *construct_candidate(candidate_type(edge));
	
	model.apply_predict(candidate.state, node_states, candidate.out_edge, candidate.out_edge.features, true);

	//increment_attribute(candidate.out_edge.attributes, attr_predict);

	candidate.score = function(candidate.out_edge.features);
	
	//int cardinality = cicada::semiring::log(counts_inside[graph.goal]) - cicada::semiring::log(counts_outside[edge.head]);
	int cardinality = cicada::semiring::log(counts_inside[graph.goal]);
	
	edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	for (edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	  cardinality -= cicada::semiring::log(counts_inside[*titer]);
	
	// - 1 to make an adjustment...
	buf[cardinality - 1].push(&candidate);
      }
    }
    
    void process_bins(const int step, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      //std::cerr << "step: " << step << " buf: " << buf[step].size() << std::endl;

      predictions.clear();

      for (size_type num_pop = 0; ! buf[step].empty() && num_pop != pop_size_max; ++ num_pop) {
	const candidate_type* item_top = buf[step].top();
	candidate_type* item = const_cast<candidate_type*>(item_top);
	
	buf[step].pop();
	
#if 0
	std::cerr << "popped head: " << item->in_edge->head
		  << " edge: " << item->in_edge->id
		  << " rule: " << *(item->in_edge->rule)
		  << " dot: " << item->dot
		  << std::endl;
#endif
	
	// we will iterate until completion...
	for (int iter = 0;; ++ iter) {
#if 0
	  std::cerr << "\titer: " << iter
		    << " head: " << item->in_edge->head
		    << " edge: " << item->in_edge->id
		    << " rule: " << *(item->in_edge->rule)
		    << " dot: " << item->dot
		    << std::endl;
#endif	    
	  
	  const bool is_goal = (graph_in.goal == item->in_edge->head);
	  
 	  const rule_type::symbol_set_type& target = item->in_edge->rule->rhs;
	  
	  // scan... and score... state will be updated...
	  if (! target.empty() && item->dot < static_cast<int>(target.size()) && ! target[item->dot].is_non_terminal()) {
	    //std::cerr << "\t\tscan" << std::endl;

	    // new item...
	    if (item == item_top) {
	      item = construct_candidate(*item);
	      item->state = model.clone(item->state);
	    }

	    feature_set_type features;
	    
	    model.apply_scan(item->state, node_states, item->out_edge, item->dot, features, is_goal);

	    //increment_attribute(item->out_edge.attributes, attr_scan);
	    
	    item->out_edge.features += features;
	    item->score *= function(features);
	    
	    // proceed dot(s)
	    for (/**/; item->dot < static_cast<int>(target.size()) && ! target[item->dot].is_non_terminal(); ++ item->dot);
	  }
	  
	  if (item->dot == static_cast<int>(target.size())) {
	    // complete...
	    //std::cerr << "\t\tcompletion" << std::endl;
	    
	    // new item...
	    if (item == item_top) {
	      item = construct_candidate(*item);
	      item->state = model.clone(item->state);
	    }
	    
	    feature_set_type features;
	    
	    model.apply_complete(item->state, node_states, item->out_edge, features, is_goal);

	    //increment_attribute(item->out_edge.attributes, attr_complete);
	    
	    item->out_edge.features += features;
	    item->score *= function(features);
	    
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
	      destroy_candidate(item);
	      
	      break;
	    } else {
	      // we will merge by the parent and state
	      std::pair<typename state_node_map_type::iterator, bool> result = states[item->in_edge->head].insert(std::make_pair(std::make_pair(item->parent, item->state), 0));
	      if (result.second) {
		node_states.push_back(item->state);
		result.first->second = graph_out.add_node().id;
	      } else
		model.deallocate(item->state);
	      
	      edge_type& edge_new = graph_out.add_edge(item->out_edge);
	      graph_out.connect_edge(edge_new.id, result.first->second);
	      
	      // we will not propagate...
	      if (! result.second) {
		// update score in the parent!
		// I is possible since parent is alwas popped from the priority queues... is it correct?
		const_cast<score_type&>(item->parent->score) = std::max(item->parent->score, item->score);
		
		destroy_candidate(item);
		
		break;
	      }
	      
	      // some trick:
	      // item->state is either deleted or inserted in states[item->in_edge->head].nodes
	      // thus, we simply copy stat from item->parent...
	      // but reassigned from siter->first by cloning...
	      
	      const state_type state = model.clone(result.first->first.second);
	      const score_type score = item->score;
	      
	      // we copy from parent, and use the score/state from the current item
	      *item = *(item->parent);
	      item->state = state;
	      item->score = score;
	      
	      const rule_type::symbol_set_type& target = item->in_edge->rule->rhs;
	      
	      const int __non_terminal_index = target[item->dot].non_terminal_index();
	      const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, item->dot_antecedent, __non_terminal_index - 1);
	      
	      item->out_edge.tails[antecedent_index] = result.first->second;
	      
	      ++ item->dot;
	      ++ item->dot_antecedent;
	      
	      continue;
	    }
	    
	  } else {
	    //std::cerr << "\t\tprediction" << std::endl;

	    // we will uniquify predictions...
	    std::pair<typename candidate_unique_set_type::iterator, bool> result = predictions.insert(item);
	    if (! result.second) {
#if 0
	      std::cerr << "multiple prediction head: " << item->in_edge->head
			<< " edge: " << item->in_edge->id
			<< " rule: " << *(item->in_edge->rule)
			<< " dot: " << item->dot
			<< std::endl;
#endif
	      
#if 0
	      if (! model_type::state_equal(model.state_size())((*result.first)->state, item->state))
		std::cerr << "merged but different state?" << std::endl;
#endif
	      
	      // re-assign score
	      const_cast<score_type&>((*result.first)->score) = std::max((*result.first)->score, item->score);
	      
	      // this model is no longer used...
	      model.deallocate(item->state);
	      
	      // destroy this candidate
	      destroy_candidate(item);
	    }
	    
	    break;
	  }
	}
      }

      //std::cerr << "predicted: " << predictions.size() << std::endl;
      
      // perform actual prediction!
      typename candidate_unique_set_type::const_iterator piter_end = predictions.end();
      for (typename candidate_unique_set_type::const_iterator piter = predictions.begin(); piter != piter_end; ++ piter) {
	const candidate_type& parent = *(*piter);
	const rule_type::symbol_set_type& target = parent.in_edge->rule->rhs;
	
	const int __non_terminal_index = target[parent.dot].non_terminal_index();
	const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, parent.dot_antecedent, __non_terminal_index - 1);
	
	const node_type& antecedent_node = graph_in.nodes[parent.out_edge.tails[antecedent_index]];
		
	node_type::edge_set_type::const_iterator eiter_end = antecedent_node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = antecedent_node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const edge_type& edge = graph_in.edges[*eiter];
	  
	  candidate_type& candidate = *construct_candidate(candidate_type(parent, edge));
	  
	  candidate.state = model.clone(parent.state);
	  
	  model.apply_predict(candidate.state, node_states, candidate.out_edge, candidate.out_edge.features, false);
	  
	  //increment_attribute(candidate.out_edge.attributes, attr_predict);
	  
	  candidate.score = parent.score * function(candidate.out_edge.features);
	  
	  //int cardinality = cicada::semiring::log(counts_inside[graph_in.goal]) - cicada::semiring::log(counts_outside[edge.head]);
	  int cardinality = cicada::semiring::log(counts_inside[edge.head]);
	  
	  edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	  for (edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	    cardinality -= cicada::semiring::log(counts_inside[*titer]);
	  
	  const size_type step_next = step + cardinality;
	  
#if 0
	  // this is not necessary...
	  if (step_next >= buf.size())
	    buf.resize(step_next + 1);
#endif
	  
	  buf[step_next].push(&candidate);
	}
      }
      
      // clear buf...
      // since we will not propagate the states in the buf, we will deallocate...
      // However, we should not reclaim candidate since they may be pointed out by other predicted candidates
      typename candidate_heap_type::const_iterator biter_end = buf[step].end();
      for (typename candidate_heap_type::const_iterator biter = buf[step].begin(); biter != biter_end; ++ biter)
	model.deallocate((*biter)->state);
      
      // shrink wrap...
      buf[step].clear();
      candidate_heap_type(buf[step]).swap(buf[step]);
    };
    
    candidate_type* construct_candidate(const candidate_type& cand)
    {
      if (candidate_list) {
	candidate_type* cand_cached = candidate_list;
	candidate_list = const_cast<candidate_type*>(cand_cached->parent);
	
	// assignment!
	*cand_cached = cand;
	
	return cand_cached;
      }
      
      candidates.push_back(cand);
      return &candidates.back();
    }

    void destroy_candidate(const candidate_type* cand)
    {
      if (! cand) return;
      
      const_cast<candidate_type*>(cand)->parent = candidate_list;
      candidate_list = const_cast<candidate_type*>(cand);
    }

#if 0
    void increment_attribute(attribute_set_type& attributes, const attribute_type& attr)
    {
      attribute_set_type::iterator iter = attributes.find(attr);
      if (iter == attributes.end())
	attributes[attr] = attribute_set_type::int_type(1);
      else
	attributes[attr] = boost::apply_visitor(__attribute_integer(), iter->second) + 1;
    }
#endif
    
  private:
    candidate_set_type  candidates;
    candidate_type*     candidate_list;
    state_set_type      node_states;
    cand_state_set_type states;

    candidate_heap_set_type buf;
    
    count_set_type counts_inside;
    //count_set_type counts_outside;
    candidate_unique_set_type predictions;

    const model_type& model;
    const function_type& function;
    size_type  pop_size_max;

#if 0
    attribute_type attr_scan;
    attribute_type attr_complete;
    attribute_type attr_predict;
#endif
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
