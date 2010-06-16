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

  private:
    typedef std::vector<feature_function_ptr_type, std::allocator<feature_function_ptr_type > > model_set_type;

  public:
    typedef model_set_type::reference       reference;
    typedef model_set_type::const_reference const_reference;

    typedef model_set_type::const_iterator iterator;
    typedef model_set_type::const_iterator const_iterator;
    
  public:
    Model() {}

  public:
    
    state_type operator()(const hypergraph_type& graph,
			  const state_set_type& node_states,
			  edge_type& edge,
			  feature_set_type& estimates) const;
    
    void operator()(const state_type& state,
		    edge_type& edge) const;
    
    const_reference operator[](size_type pos) const { return models[pos]; }
    reference operator[](size_type pos) { return models[pos]; }

    const_iterator begin() const { return models.begin(); }
    iterator begin() { return models.begin(); }

    const_iterator end() const { return models.end(); }
    iterator end() { return models.end(); }
    
    // you should call this at least once, when you are going to use this model.
    void initialize();

    bool empty() const { return models.empty(); }
    size_type size() const { return models.size(); }
    
    void push_back(const feature_function_ptr_type& x)
    {
      models.push_back(x);
    }
    
    void clear()
    {
      models.clear();
      offsets.clear();
      states_size = 0;
    }

    bool is_stateless() const { return states_size == 0; }
    
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
