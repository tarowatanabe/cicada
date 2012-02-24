// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__APPLY_EXACT__HPP__
#define __CICADA__APPLY_EXACT__HPP__ 1

#include <vector>

#include <cicada/apply_state_less.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/model.hpp>
#include <cicada/semiring/traits.hpp>

#include <utils/dense_hash_map.hpp>
#include <utils/small_vector.hpp>
#include <utils/hashmurmur.hpp>

namespace cicada
{
  
  // a naive algorithm...


  struct ApplyExact
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
        
    typedef utils::small_vector<int, std::allocator<int> > index_set_type;
    
    typedef std::vector<id_type, std::allocator<id_type> > node_set_type;
    typedef std::vector<node_set_type, std::allocator<node_set_type> > node_map_type;
    
    typedef google::dense_hash_map<state_type, id_type, model_type::state_hash, model_type::state_equal > state_node_map_type;
        
    ApplyExact(const model_type& _model)
      : model(_model)
    {  }
    
    void operator()(const hypergraph_type& graph_in,
		    hypergraph_type&       graph_out)
    {
      const_cast<model_type&>(model).initialize();

      if (model.is_stateless()) {
	ApplyStateLess __applier(model);
	__applier(graph_in, graph_out);
      } else {
	node_map.clear();
	node_map.reserve(graph_in.nodes.size());
	node_map.resize(graph_in.nodes.size());
	
	node_states.clear();
	node_states.reserve(graph_in.nodes.size() * 10000);
	
	graph_out.clear();
	for (id_type node_id = 0; node_id < graph_in.nodes.size(); ++ node_id)
	  process(node_id, graph_in, graph_out);
	
	// topologically sort...
	graph_out.topologically_sort();
	
	// re-initialize again...
	const_cast<model_type&>(model).initialize();
      }
    };
    
  private:
    
    void process(id_type v, const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      const node_type& node = graph_in.nodes[v];
      const bool is_goal(v == graph_in.goal);

      state_node_map_type buf(node.edges.size(), model_type::state_hash(model.state_size()), model_type::state_equal(model.state_size()));
      buf.set_empty_key(state_type());
      
      node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = graph_in.edges[*eiter];
	
	index_set_type j_ends(edge.tails.size(), 0);
	index_set_type j(edge.tails.size(), 0);

	for (size_t i = 0; i != edge.tails.size(); ++ i)
	  j_ends[i] = node_map[edge.tails[i]].size();

	edge_type::node_set_type tails(edge.tails.size());
	
	for (;;) {
	  
	  for (size_t i = 0; i != edge.tails.size(); ++ i)
	    tails[i] = node_map[edge.tails[i]][j[i]];
	  
	  edge_type& edge_new = graph_out.add_edge(tails.begin(), tails.end());
	  edge_new.head = v;
	  edge_new.rule = edge.rule;
	  edge_new.features   = edge.features;
	  edge_new.attributes = edge.attributes;

	  const state_type state = model.apply(node_states, edge_new, edge_new.features, is_goal);
	  
	  // hypothesis recombination
	  
	  if (is_goal) {
	    if (graph_out.goal == hypergraph_type::invalid) {
	      graph_out.goal = graph_out.add_node().id;
	      node_states.push_back(state);
	    } else
	      model.deallocate(state);
	    
	    node_type& node = graph_out.nodes[graph_out.goal];
	    
	    graph_out.connect_edge(edge_new.id, node.id);
	  } else {
	    typedef std::pair<state_node_map_type::iterator, bool > result_type;
	    
	    result_type result = buf.insert(std::make_pair(state, 0));
	    if (result.second) {
	      result.first->second = graph_out.add_node().id;
	      
	      node_states.push_back(state);
	      
	      node_map[edge.head].push_back(result.first->second);
	    } else
	      model.deallocate(state);
	    
	    graph_out.connect_edge(edge_new.id, result.first->second);
	  }
	  
	  // proceed to the next id...

	  size_t index = 0;
	  for (/**/; index != edge.tails.size(); ++ index) {
	    ++ j[index];
	    if (j[index] < j_ends[index]) break;
	    j[index] = 0;
	  }
	  
	  // finished!
	  if (index == edge.tails.size()) break;
	}
      }
    }
    
  private:
    node_map_type       node_map;
    state_set_type      node_states;
    
    const model_type& model;
  };


  inline
  void apply_exact(const Model& model, const HyperGraph& source, HyperGraph& target)
  {
    ApplyExact __apply(model);

    __apply(source, target);
  }
  
  inline
  void apply_exact(const Model& model, HyperGraph& source)
  {
    HyperGraph target;
    
    ApplyExact __apply(model);
    
    __apply(source, target);
    
    source.swap(target);
  }

};

#endif
