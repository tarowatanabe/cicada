// -*- mode: c++ -*-

// a model composed of many feature functions


#ifndef __CICADA__MODEL__HPP__
#define __CICADA__MODEL__HPP__ 1

#include <cicada/feature_function.hpp>
#include <cicada/hypergraph.hpp>

#include <boost/shared_ptr.hpp>

#include <utils/symbol.hpp>

namespace cicada
{
  class Model
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef float     weight_type;

    typedef cicada::Symbol     symbol_type;
    typedef cicada::Feature    feature_type;
    typedef cicada::HyperGraph hypergraph_type;
    
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    typedef hypergraph_type::feature_set_type feature_set_type;
    
  public:
    typedef utils::simple_vector<char, std::allocator<char> > state_type;
    
    typedef std::vector<state_type, std::allocator<state_type> > state_set_type;
    
    typedef FeatureFunction                          feature_function_type;
    typedef boost::shared_ptr<feature_function_type> feature_function_ptr_type;
    typedef std::vector<feature_function_ptr_type, std::allocator<feature_function_ptr_type > > model_set_type;
    
  public:
    Model() : states_size(0) {}
    Model(const model_set_type& __models)
      : models(__models),
	offsets(__models.size()),
	states_size(0) { initialize(); }

  public:
    
    state_type operator()(const hypergraph_type& graph,
			  const state_set_type& node_states,
			  edge_type& edge,
			  feature_set_type& estimates) const;

    void operator()(const state_type& state,
		    edge_type& edge) const;

  private:
    void initialize();
    
  private:
    typedef std::vector<size_type, std::allocator<size_type> > offset_set_type;
    
  private:
    model_set_type  models;
    
    // offsets used to compute states...
    offset_set_type offsets;
    int             states_size;
  };
};


#endif
