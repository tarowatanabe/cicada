// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__APPLY_CUBE_PRUNE__HPP__
#define __CICADA__APPLY_CUBE_PRUNE__HPP__ 1

#include <cicada/apply_state_less.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/model.hpp>

#include <cicada/semiring/traits.hpp>

#include <google/dense_hash_set>
#include <google/dense_hash_map>

#include <utils/simple_vector.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>

#include <utils/b_heap.hpp>
#include <utils/std_heap.hpp>

namespace cicada
{
  
  // implementation of 
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
  struct ApplyCubePrune
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
      const edge_type* in_edge;
      edge_type        out_edge;
      
      state_type state;
      
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
    typedef std::vector<node_score_list_type, std::allocator<node_score_list_type> > node_score_set_type;
    
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
    
    struct compare_estimate_type
    {
      // we will use greater, so that simple sort will yield estimated score order...
      bool operator()(const candidate_type* x, const candidate_type* y) const
      {
	return x->estimate > y->estimate;
      }

      bool operator()(const node_score_type& x, const node_score_type& y) const
      {
	return x.estimate > y.estimate;
      }
    };

    //typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_type;
    
    typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_base_type;
    //typedef utils::b_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type, 512 / sizeof(const candidate_type*)> candidate_heap_type;
    typedef utils::std_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type> candidate_heap_type;
    
    typedef google::dense_hash_map<state_type, candidate_type*, model_type::state_hash, model_type::state_equal > state_node_map_type;
    typedef google::dense_hash_set<const candidate_type*, candidate_hash_type, candidate_equal_type > candidate_set_unique_type;
    
    ApplyCubePrune(const model_type& _model,
		   const function_type& _function,
		   const int _cube_size_max)
      : model(_model),
	function(_function),
	cube_size_max(_cube_size_max)
    { 
      cand_unique.set_empty_key(0);
    }
    
    void operator()(const hypergraph_type& graph_in,
		    hypergraph_type&       graph_out)
    {
      const_cast<model_type&>(model).initialize();

      if (model.is_stateless()) {
	ApplyStateLess __applier(model);
	__applier(graph_in, graph_out);
      } else {
	candidates.clear();
	
	D.clear();
	D.reserve(graph_in.nodes.size());
	D.resize(graph_in.nodes.size());
	
	node_states.clear();
	node_states.reserve(graph_in.nodes.size() * cube_size_max);
	
	graph_out.clear();
	for (id_type node_id = 0; node_id < graph_in.nodes.size(); ++ node_id)
	  kbest(node_id, graph_in, graph_out);
      
	// topologically sort...
	graph_out.topologically_sort();
	
	// re-initialize again...
	const_cast<model_type&>(model).initialize();
      }
    };
    
  private:
    
    void kbest(id_type v, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      //std::cerr << "kbest node: " << v << std::endl;
      
      // clear candidates!
      candidates.clear();
      
      const node_type& node = graph_in.nodes[v];
      const bool is_goal(v == graph_in.goal);
      
      // for each incoming e, cand \leftarrow { <e, 1>}
      
      cand.clear();
      cand.reserve(node.edges.size() * cube_size_max);
      
      cand_unique.clear();
      
      node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = graph_in.edges[*eiter];
	const index_set_type j(edge.tails.size(), 0);
	
	const candidate_type* item = make_candidate(edge, j, graph_out, is_goal);
	
	cand.push(item);
	cand_unique.insert(item);
      }
      
      //std::cerr << "perform cube-prune" << std::endl;
      
      state_node_map_type buf(cand.size(), model_type::state_hash(model.state_size()), model_type::state_equal(model.state_size()));
      buf.set_empty_key(state_type());
      
      for (size_type num_pop = 0; !cand.empty() && num_pop != cube_size_max; ++ num_pop) {
	// pop-best...
	const candidate_type* item = cand.top();
	cand.pop();
	
	push_succ(*item, is_goal, cand, cand_unique, graph_out);
	append_item(*item, is_goal, buf, graph_out);
      }
      
      //std::cerr << "finished" << std::endl;
      
      // sort buf to D(v)
      D[v].reserve(buf.size());
      
      typename state_node_map_type::const_iterator biter_end = buf.end();
      for (typename state_node_map_type::const_iterator biter = buf.begin(); biter != biter_end; ++ biter)
	D[v].push_back(node_score_type(biter->second->out_edge.head, biter->second->score, biter->second->estimate));
      
      std::sort(D[v].begin(), D[v].end(), compare_estimate_type());

      typename candidate_heap_type::const_iterator hiter_end = cand.end();
      for (typename candidate_heap_type::const_iterator hiter = cand.begin(); hiter != hiter_end; ++ hiter)
	model.deallocate((*hiter)->state);
    }
    
    // append item to buf
    void append_item(const candidate_type& item,
		     const bool is_goal,
		     state_node_map_type& buf,
		     hypergraph_type& graph)
    {
      edge_type& edge_new = graph.add_edge(item.out_edge);

#if 0
      std::cerr << "edge-id: " << edge_new.id
		<< " head: " << edge_new.head
		<< std::endl;
#endif
      
      if (is_goal) {
	// perform hypothesis re-combination toward goal-node...
	
	if (graph.goal == hypergraph_type::invalid) {
	  graph.goal = graph.add_node().id;
	  node_states.push_back(item.state);
	} else
	  model.deallocate(item.state);
	
	node_type& node = graph.nodes[graph.goal];
	
	graph.connect_edge(edge_new.id, node.id);
      } else {
	// hypothesis re-combination...
	
	typedef std::pair<typename state_node_map_type::iterator, bool> result_type;
	
	result_type result = buf.insert(std::make_pair(item.state, const_cast<candidate_type*>(&item)));
	
	if (result.second) {
	  //std::cerr << "added node!" << std::endl;
	  
	  result.first->second->out_edge.head = graph.add_node().id;
	  node_states.push_back(item.state);
	} else
	  model.deallocate(item.state);
	
	//std::cerr << "node: " << biter->second->node << std::endl;
	
	candidate_type& item_graph = *(result.first->second);
	
	node_type& node = graph.nodes[item_graph.out_edge.head];
	
	graph.connect_edge(edge_new.id, node.id);
	
	// check if we found better derivation.. 
	// it may happen due to insufficient contextual information during parsing...
	if (item.score > item_graph.score) {
	  item_graph.score  = item.score;
	  item_graph.estimate = item.estimate;
	}
      }
      
      //std::cerr << "finished append-item" << std::endl;
    }
    
    // push succ...
    size_type push_succ(const candidate_type& candidate,
			const bool is_goal,
			candidate_heap_type& cand,
			candidate_set_unique_type& candidates_unique,
			hypergraph_type& graph_out)
    {
      // proceed to the next until we find better candidate...
      // @InProceedings{iglesias-EtAl:2009:EACL,
      //  author    = {Iglesias, Gonzalo  and  de Gispert, Adri\`{a}  and  Banga, Eduardo R.  and  Byrne, William},
      //  title     = {Rule Filtering by Pattern for Efficient Hierarchical Translation},
      //  booktitle = {Proceedings of the 12th Conference of the European Chapter of the ACL (EACL 2009)},
      //  month     = {March},
      //  year      = {2009},
      //  address   = {Athens, Greece},
      //  publisher = {Association for Computational Linguistics},
      //  pages     = {380--388},
      //  url       = {http://www.aclweb.org/anthology/E09-1044}
      // }

      candidate_type query(candidate.j);
      
      query.in_edge = candidate.in_edge;
      index_set_type& j = query.j;
      
      size_type inserted = 0;
      for (size_t i = 0; i != candidate.j.size(); ++ i) {
	
	const int j_i_prev = j[i];
	++ j[i];
	
	for (/**/; j[i] < static_cast<int>(D[candidate.in_edge->tails[i]].size()); ++ j[i]) {
	  
	  if (candidates_unique.find(&query) == candidates_unique.end()) {
	    // new candidate...
	    const candidate_type* candidate_new = make_candidate(*candidate.in_edge, j, graph_out, is_goal);
	    
	    cand.push(candidate_new);
	    candidates_unique.insert(candidate_new);
	    
	    ++ inserted;
	    
	    break;
	  }
	}
	
	j[i] = j_i_prev;
      }
      
      //std::cerr << "inserted: " << inserted << std::endl;
      
      return inserted;
    }
    
    const candidate_type* make_candidate(const edge_type& edge, const index_set_type& j, const hypergraph_type& graph, const bool is_goal)
    {
      //std::cerr << "make candidate for: " << *(edge.rule) << std::endl;

      candidates.push_back(candidate_type(edge, j));
      
      candidate_type& candidate = candidates.back();
      
      candidate.out_edge.tails = edge_type::node_set_type(j.size());
      
      candidate.score = semiring::traits<score_type>::one();
      candidate.estimate = semiring::traits<score_type>::one();
      for (size_t i = 0; i != j.size(); ++ i) {
	const node_score_type& antecedent = D[edge.tails[i]][j[i]];
	
	candidate.out_edge.tails[i] = antecedent.node;
	candidate.score *= antecedent.score;
      }
      
      // perform actual model application...
      
      feature_set_type estimates;
      candidate.state = model.apply(node_states, candidate.out_edge, candidate.out_edge.features, estimates, is_goal);
      
      candidate.score    *= function(candidate.out_edge.features);
      candidate.estimate *= function(estimates);
      candidate.estimate *= candidate.score;

      //std::cerr << "make candidate done" << std::endl;
      
      return &candidate;
    };
    
  private:
    candidate_set_type  candidates;
    node_score_set_type D;
    state_set_type      node_states;

    candidate_heap_type       cand;
    candidate_set_unique_type cand_unique;

    const model_type& model;
    const function_type& function;
    size_type  cube_size_max;
  };
  
  template <typename Function>
  inline
  void apply_cube_prune(const Model& model, const HyperGraph& source, HyperGraph& target, const Function& func, const int cube_size)
  {
    ApplyCubePrune<typename Function::value_type, Function>(model, func, cube_size)(source, target);
  }

  template <typename Function>
  inline
  void apply_cube_prune(const Model& model, HyperGraph& source, const Function& func, const int cube_size)
  {
    HyperGraph target;
    
    ApplyCubePrune<typename Function::value_type, Function>(model, func, cube_size)(source, target);
    
    source.swap(target);
  }

};

#endif
