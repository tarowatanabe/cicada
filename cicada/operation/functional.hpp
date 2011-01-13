// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
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
    
    struct shortest_length_function
    {
      typedef cicada::Vocab vocab_type;
      typedef cicada::Rule rule_type;
      typedef cicada::semiring::Tropical<int> value_type;
	
      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	int length = 0;
	rule_type::symbol_set_type::const_iterator siter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator siter = edge.rule->rhs.begin(); siter != siter_end; ++ siter)
	  length += (*siter != vocab_type::EPSILON && siter->is_terminal());
	  
	// since we will "max" at operator+, we will collect negative length
	return cicada::semiring::traits<value_type>::log(- length);
      }
    };

    struct longest_length_function
    {
      typedef cicada::Vocab vocab_type;
      typedef cicada::Rule rule_type;
      typedef cicada::semiring::Tropical<int> value_type;
	
      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	int length = 0;
	rule_type::symbol_set_type::const_iterator siter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator siter = edge.rule->rhs.begin(); siter != siter_end; ++ siter)
	  length += (*siter != vocab_type::EPSILON && siter->is_terminal());
	  
	// since we will "max" at operator+, we will collect positive length
	return cicada::semiring::traits<value_type>::log(length);
      }
    };

  };
};

#endif
