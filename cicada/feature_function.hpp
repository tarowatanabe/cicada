// -*- mode: c++ -*-

#ifndef __CICADA__FEATURE_FUNCTION__HPP__
#define __CICADA__FEATURE_FUNCTION__HPP__ 1

// we will define a basic feature functions
// stateless feature will have zero size feature

#include <vector>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/feature_vector.hpp>
#include <cicada/feature.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  class FeatureFunction
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef cicada::Symbol     symbol_type;
    typedef cicada::Vocab      vocab_type;
    typedef cicada::Feature    feature_type;
    typedef cicada::HyperGraph hypergraph_type;
    typedef cicada::Rule       rule_type;
    
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    typedef hypergraph_type::feature_set_type feature_set_type;
    
  public:
    typedef void* state_ptr_type;
    typedef std::vector<state_ptr_type, std::allocator<state_ptr_type> > state_ptr_set_type;

    typedef FeatureFunction                          feature_function_type;
    typedef boost::shared_ptr<feature_function_type> feature_function_ptr_type;
    
  public:
    FeatureFunction()
      : __state_size(0) {}
    FeatureFunction(int state_size)
      : __state_size(state_size) {}
    FeatureFunction(int state_size, const feature_type& feature_name)
      : __state_size(state_size), __feature_name(feature_name) {}
    virtual ~FeatureFunction() {}

  public:
    
    static feature_function_ptr_type create(const std::string& parameter);
    static std::string               lists();
    
  public:
    virtual void operator()(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features,
			    feature_set_type& estimates) const {}
    virtual void operator()(const state_ptr_type& state,
			    feature_set_type& features) const {}
    
        
    size_type state_size() const { return __state_size; }
    const feature_type& feature_name() const { return __feature_name; }

  protected:
    size_type    __state_size;
    feature_type __feature_name;
  };
  
};

#endif
