// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__FUNCTIONAL__HPP__
#define __CICADA__OPERATION__FUNCTIONAL__HPP__ 1

#include <cicada/operation.hpp>
#include <cicada/semiring.hpp>

namespace cicada
{
  namespace operation
  {
    template <typename Weight>
    struct weight_scaled_function
    {
      typedef cicada::Operation::hypergraph_type hypergraph_type;
      typedef cicada::Operation::weight_set_type weight_set_type;
      typedef Weight value_type;
  
      weight_scaled_function(const weight_set_type& __weights, const double& __scale)
	: weights(__weights), scale(__scale) {}
      
      const weight_set_type& weights;
      const double scale;
  
      value_type operator()(const hypergraph_type::edge_type& x) const
      {
	return cicada::semiring::traits<value_type>::log(x.features.dot(weights) * scale);
      }
      
      template <typename FeatureSet>
      value_type operator()(const FeatureSet& x) const
      {
	return cicada::semiring::traits<value_type>::log(x.dot(weights) * scale);
      }
    };


    template <typename Weight>
    struct weight_scaled_function_one
    {
      typedef cicada::Operation::hypergraph_type hypergraph_type;
      typedef cicada::Operation::weight_set_type weight_set_type;
      typedef Weight value_type;
  
      weight_scaled_function_one(const double& __scale)
	: scale(__scale) {}
      
      const double scale;
      
      value_type operator()(const hypergraph_type::edge_type& x) const
      {
	return cicada::semiring::traits<value_type>::log(x.features.sum() * scale);
      }
      
      template <typename FeatureSet>
      value_type operator()(const FeatureSet& x) const
      {
	return cicada::semiring::traits<value_type>::log(x.sum() * scale);
      }
    };

    template <typename Weight>
    struct weight_function
    {
      typedef cicada::Operation::hypergraph_type hypergraph_type;
      typedef cicada::Operation::weight_set_type weight_set_type;
      typedef Weight value_type;

      weight_function(const weight_set_type& __weights)
	: weights(__weights) {}

      const weight_set_type& weights;

      value_type operator()(const hypergraph_type::edge_type& x) const
      {
	return cicada::semiring::traits<value_type>::log(x.features.dot(weights));
      }
  
      template <typename FeatureSet>
      value_type operator()(const FeatureSet& x) const
      {
	return cicada::semiring::traits<value_type>::log(x.dot(weights));
      }
    };

    template <typename Weight>
    struct weight_function_one
    {
      typedef cicada::Operation::hypergraph_type hypergraph_type;
      typedef cicada::Operation::weight_set_type weight_set_type;
      typedef Weight value_type;
      
      value_type operator()(const hypergraph_type::edge_type& x) const
      {
	return cicada::semiring::traits<value_type>::log(x.features.sum());
      }
  
      template <typename FeatureSet>
      value_type operator()(const FeatureSet& x) const
      {
	return cicada::semiring::traits<value_type>::log(x.sum());
      }
    };
    
  };
};

#endif
