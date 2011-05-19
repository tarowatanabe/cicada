//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_IMPL__HPP__
#define __CICADA_LEARN_IMPL__HPP__ 1

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"
#include "cicada/hypergraph.hpp"

#include "cicada/inside_outside.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/dot_product.hpp"
#include "cicada/semiring.hpp"

typedef boost::filesystem::path path_type;

typedef cicada::Symbol          symbol_type;
typedef cicada::Vocab           vocab_type;
typedef cicada::Sentence        sentence_type;
typedef cicada::Lattice         lattice_type;
typedef cicada::Rule            rule_type;
typedef cicada::HyperGraph      hypergraph_type;

typedef hypergraph_type::feature_set_type    feature_set_type;
typedef feature_set_type::feature_type feature_type;
typedef cicada::WeightVector<double>   weight_set_type;

struct OptimizerBase
{
  OptimizerBase() {}
  
  typedef cicada::semiring::Log<double> weight_type;
  typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
  typedef cicada::WeightVector<weight_type, std::allocator<weight_type> >  expectation_type;
};

struct OptimizerSGDL2 : public OptimizerBase
{
  OptimizerSGDL2(const size_t& __instances,
		 const double& C) : instances(__instances), samples(0), epoch(0), lambda(C / __instances), weight_scale(1.0), weight_norm(0.0) {}
  
  void initialize()
  {
    samples = 0;
    
    weight_scale = 1.0;
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
    
    objective = 0.0;
  }
  void finalize()
  {
    weights *= weight_scale;
    
    weight_scale = 1.0;
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
  }
  
  void operator()(const gradient_type& correct, 
		  const gradient_type& gradient, 
		  const weight_type& Z_correct,
		  const weight_type& Z)
  {
    //const double eta = 1.0 / (1.0 + double(epoch) / graphs.size());
    //const double eta = 1.0 / (lambda * (epoch + 2));
    const double eta = 0.2 * std::pow(0.85, double(epoch) / instances);
    ++ epoch;

    rescale(1.0 - eta * lambda);
    
    gradient_type::const_iterator citer_end = correct.end();
    for (gradient_type::const_iterator citer = correct.begin(); citer != citer_end; ++ citer)
      update(weights[citer->first], citer->second * eta);
    
    gradient_type::const_iterator miter_end = gradient.end();
    for (gradient_type::const_iterator miter = gradient.begin(); miter != miter_end; ++ miter)
      update(weights[miter->first], - miter->second * eta);

    // projection...
    if (weight_norm > 1.0 / lambda)
      rescale(std::sqrt(1.0 / (lambda * weight_norm)));
    
    objective += double(log(Z_correct) - log(Z)) * weight_scale;
    
    if (weight_scale < 0.01 || 100 < weight_scale) {
      weights *= weight_scale;
      weight_scale = 1.0;
      weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
    }
    
    ++ samples;
  }

  void update(double& x, const double& alpha)
  {
    weight_norm += 2.0 * x * alpha * weight_scale + alpha * alpha;
    x += alpha / weight_scale;
  }
  
  void rescale(const double scaling)
  {
    weight_norm *= scaling * scaling;
    if (scaling != 0.0)
      weight_scale *= scaling;
    else {
      weight_scale = 1.0;
      std::fill(weights.begin(), weights.end(), 0.0);
    }
  }
  
  size_t instances;
  size_t samples;
  size_t epoch;
  double lambda;

  double weight_scale;
  double weight_norm;
  
  double objective;
  weight_set_type weights;
};


struct OptimizerSGDL1 : public OptimizerBase
{
  typedef cicada::WeightVector<double> penalty_set_type;

  OptimizerSGDL1(const size_t& __instances, const double& C)
    : instances(__instances), samples(0), epoch(0), lambda(C / __instances), penalties(), penalty(0.0), weight_scale(1.0) {}

  void initialize()
  {
    samples = 0;
    objective = 0.0;
    weight_scale = 1.0;
  }
  void finalize()
  {

  }
  
  void operator()(const gradient_type& correct, 
		  const gradient_type& gradient, 
		  const weight_type& Z_correct,
		  const weight_type& Z)
  {
    //const double eta = 1.0 / (1.0 + double(epoch) / graphs.size());
    //const double eta = 1.0 / (lambda * (epoch + 2));
    const double eta = 0.2 * std::pow(0.85, double(epoch) / instances);
    ++ epoch;
    
    gradient_type::const_iterator citer_end = correct.end();
    for (gradient_type::const_iterator citer = correct.begin(); citer != citer_end; ++ citer) {
      weights[citer->first] += eta * citer->second;
      
      apply(weights[citer->first], penalties[citer->first], penalty);
    }
    
    gradient_type::const_iterator miter_end = gradient.end();
    for (gradient_type::const_iterator miter = gradient.begin(); miter != miter_end; ++ miter) {
      weights[miter->first] -= eta * miter->second;
      
      apply(weights[miter->first], penalties[miter->first], penalty);
    }
    
    objective += double(log(Z_correct) - log(Z));
    ++ samples;
  }
  
  void apply(double& x, double& penalty, const double& cummulative)
  {
    const double x_half = x;
    if (x > 0.0)
      x = std::max(0.0, x - penalty - cummulative);
    else if (x < 0.0)
      x = std::min(0.0, x - penalty + cummulative);
    penalty += x - x_half;
  }
  
  size_t instances;
  size_t samples;
  size_t epoch;
  double lambda;
  
  penalty_set_type penalties;
  double penalty;
  
  double objective;
  weight_set_type weights;
  double weight_scale;
};


#endif
