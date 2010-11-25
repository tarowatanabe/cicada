// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__APPLY_STATE_LESS__HPP__
#define __CICADA__APPLY_STATE_LESS__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/model.hpp>
#include <cicada/semiring/traits.hpp>

#include <utils/simple_vector.hpp>
#include <utils/hashmurmur.hpp>

namespace cicada
{
  
  // a naive algorithm...


  struct ApplyStateLess
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
        
    
    ApplyStateLess(const model_type& _model)
      : model(_model)
    {  }
    
    void operator()(const hypergraph_type& graph_in,
		    hypergraph_type&       graph_out)
    {
      const_cast<model_type&>(model).initialize();

      if (! model.is_stateless())
	throw std::runtime_error("we do not support state-full feature application");
      
      node_states.clear();
      node_states.reserve(graph_in.nodes.size());
      node_states.resize(graph_in.nodes.size());
      
      graph_out = graph_in;
      
      // traverse in topological order...
      
      hypergraph_type::node_set_type::const_iterator niter_end = graph_out.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph_out.nodes.begin(); niter != niter_end; ++ niter) {
	const node_type& node = *niter;
	
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  edge_type& edge = graph_out.edges[*eiter];
	  
	  feature_set_type estimates;
	  model.apply(node_states, edge, edge.features, estimates, node.id == graph_out.goal);
	}
      }
      

      // re-initialize again...
      const_cast<model_type&>(model).initialize();
    };
    
  private:
    state_set_type      node_states;
    
    const model_type& model;
  };


  inline
  void apply_state_less(const Model& model, const HyperGraph& source, HyperGraph& target)
  {
    ApplyStateLess __apply(model);

    __apply(source, target);
  }
  
  inline
  void apply_state_less(const Model& model, HyperGraph& source)
  {
    HyperGraph target;
    
    ApplyStateLess __apply(model);
    
    __apply(source, target);
    
    source.swap(target);
  }

};

#endif
