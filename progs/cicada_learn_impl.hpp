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
		 const double& C) : instances(__instances), samples(0), epoch(0), lambda(C), weight_scale(1.0), weight_norm(0.0) {}
  
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

    gradient_type::const_iterator citer = correct.begin();
    gradient_type::const_iterator citer_end = correct.end();
    
    gradient_type::const_iterator miter = gradient.begin();
    gradient_type::const_iterator miter_end = gradient.end();

    while (citer != citer_end && miter != miter_end) {
      if (citer < miter) {
	update(weights[citer->first], double(citer->second) * eta);
	++ citer;
      } else if (miter < citer) {
	update(weights[miter->first], - double(miter->second) * eta);
	++ miter;
      } else {
	const double alpha = double(citer->second) - double(miter->second);
	if (alpha != 0.0)
	  update(weights[citer->first], alpha * eta);
	++ citer;
	++ miter;
      }
    }
    
    for (/**/; citer != citer_end; ++ citer)
      update(weights[citer->first], double(citer->second) * eta);

    for (/**/; miter != miter_end; ++ miter)
      update(weights[miter->first], - double(miter->second) * eta);
    
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
    : instances(__instances), samples(0), epoch(0), lambda(C), penalties(), penalty(0.0), weight_scale(1.0) {}

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
    const double factor = 1.0 / instances;
    const double eta = 0.2 * std::pow(0.85, double(epoch) / instances);
    ++ epoch;
    
    penalty += eta * lambda;
    
    gradient_type::const_iterator citer = correct.begin();
    gradient_type::const_iterator citer_end = correct.end();
    
    gradient_type::const_iterator miter = gradient.begin();
    gradient_type::const_iterator miter_end = gradient.end();
    
    while (citer != citer_end && miter != miter_end) {
      if (citer < miter) {
	weights[citer->first] += eta * double(citer->second) * factor;
	apply(weights[citer->first], penalties[citer->first], penalty);
	
	++ citer;
      } else if (miter < citer) {
	weights[miter->first] -= eta * double(miter->second) * factor;
	apply(weights[miter->first], penalties[miter->first], penalty);
	
	++ miter;
      } else {
	weights[citer->first] += eta * (double(citer->second) - double(miter->second)) * factor;
	apply(weights[citer->first], penalties[citer->first], penalty);
	
	++ citer;
	++ miter;
      }
    }
    
    for (/**/; citer != citer_end; ++ citer) {
      weights[citer->first] += eta * double(citer->second) * factor;
      apply(weights[citer->first], penalties[citer->first], penalty);
    }
    
    for (/**/; miter != miter_end; ++ miter) {
      weights[miter->first] -= eta * double(miter->second) * factor;
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

struct OptimizerMIRA : public OptimizerBase
{
  OptimizerMIRA(const size_t& __instances,
		const double& C) : instances(__instances), samples(0), lambda(C) {}

  void initialize()
  {
    samples = 0;
    objective = 0.0;
    weight_scale = 1.0;
  }
  void finalize()
  {
    
  }
  
  template <typename Iterator>
  void operator()(Iterator first, Iterator last, const double loss)
  {
    const double margin = cicada::dot_product(weights, first, last, 0.0);
    const double variance = cicada::dot_product(first, last, first, last, 0.0);
    
    objective += loss - margin;
    
    const double alpha = std::max(0.0, std::min(1.0 / C, (loss - margin) / variance));
    
    if (alpha > 1e-10) {
      for (/**/; first != last; ++ first)
	weights[first->first] += alpha * first->second;
    
      ++ samples;
    }
  }

  void operator()(const feature_set_type& features_reward,
		  const feature_set_type& features_penalty,
		  const double loss=1.0)
  {
    const feature_set_type features(features_reward - features_penalty);
    
    const double margin = cicada::dot_product(weights, features);
    const double variance = cicada::dot_product(features, features);
    
    objective += loss - margin;
    
    const double alpha = std::max(0.0, std::min(1.0 / C, (loss - margin) / variance));
    
    if (alpha > 1e-10) {
      feature_set_type::const_iterator fiter_end = features.end();
      for (feature_set_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
	weights[fiter->first] += alpha * fiter->second;
    }
    ++ samples;
  }
  
  size_t instances;
  size_t samples;
  
  double objective;
  weight_set_type weights;
  double weight_scale;
  double lambda;
};

struct OptimizerAROW : public OptimizerBase
{
  OptimizerAROW(const size_t& __instances,
		const double& C) : instances(__instances), samples(0), lambda(C) {}

  void initialize()
  {
    samples = 0;
    objective = 0.0;
    weight_scale = 1.0;
  }
  
  void finalize()
  {
    
  }

  template <typename Iterator>
  void operator()(Iterator first, Iterator last, const double loss)
  {
    covariances.allocate(1.0);
    
    const double margin = cicada::dot_product(weights, first, last, 0.0);
    const double variance = cicada::dot_product(first, last, covariances, first, last, 0.0); // multiply covariances...
    
    objective += loss - margin;
    
    const double beta = 1.0 / (variance + C);
    const double alpha = std::max(0.0, (loss - margin) * beta);
    
    if (alpha > 1e-10) {
      for (/**/; first != last; ++ first) {
	const double var = covariances[first->first];
	
	weights[first->first]     += alpha * first->second * var;
	covariances[first->first] -= beta * (var * var) * (first->second * first->second);
      }
      
      ++ samples;
    }
  }

  void operator()(const feature_set_type& features_reward,
		  const feature_set_type& features_penalty,
		  const double loss=1.0)
  {
    feature_set_type features(features_reward - features_penalty);
    
    covariances.allocate(1.0);
    
    const double margin = cicada::dot_product(weights, features);
    const double variance = cicada::dot_product(features, covariances, features); // multiply covariances...
    
    objective += loss - margin;
    
    const double beta = 1.0 / (variance + C);
    const double alpha = std::max(0.0, (loss - margin) * beta);
    
    if (alpha > 1e-10) {
      feature_set_type::const_iterator fiter_end = features.end();
      for (feature_set_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter) {
	const double var = covariances[fiter->first];
	
	weights[fiter->first]     += alpha * fiter->second * var;
	covariances[fiter->first] -= beta * (var * var) * (fiter->second * fiter->second);
      }
    }
    ++ samples;
  }
  
  size_t instances;
  size_t samples;
  
  double objective;
  weight_set_type weights;
  double weight_scale;
  double lambda;

  weight_set_type covariances;
};

struct OptimizerCW : public OptimizerBase
{
  OptimizerCW(const size_t& __instances,
		const double& C) : instances(__instances), samples(0), lambda(C) {}

  void initialize()
  {
    samples = 0;
    objective = 0.0;
    weight_scale = 1.0;
  }
  
  void finalize()
  {
    
  }

  template <typename Iterator>
  void operator()(Iterator first, Iterator last, const double loss)
  {
    covariances.allocate(1.0);
    
    const double margin = cicada::dot_product(weights, first, last, 0.0);
    const double variance = cicada::dot_product(first, last, covariances, first, last, 0.0); // multiply covariances...
    
    objective += loss - margin;
    
    if (loss - margin > 0.0) {
      const double theta = 1.0 + 2.0 * C * (margin - loss);
      const double alpha = ((- theta + std::sqrt(theta * theta - 8.0 * C * (margin - loss - C * variance))) / (4.0 * C * variance));
      const double beta  = (2.0 * alpha * C) / (1.0 + 2.0 * alpha * C * variance);
      
      if (alpha > 1e-10 && beta > 0.0) {
	for (/**/; first != last; ++ first) {
	  const double var = covariances[first->first];
	  
	  weights[first->first]     += alpha * first->second * var;
	  covariances[first->first] -= beta * (var * var) * (first->second * first->second);
	}
	
	++ samples;
      }
    }
    
    
  }

  void operator()(const feature_set_type& features_reward,
		  const feature_set_type& features_penalty,
		  const double loss=1.0)
  {
    feature_set_type features(features_reward - features_penalty);
    
    covariances.allocate(1.0);
    
    const double margin = cicada::dot_product(weights, features);
    const double variance = cicada::dot_product(features, covariances, features); // multiply covariances...
    
    objective += loss - margin;
    
    if (loss - margin > 0.0) {
      const double theta = 1.0 + 2.0 * C * (margin - loss);
      const double alpha = ((- theta + std::sqrt(theta * theta - 8.0 * C * (margin - loss - C * variance))) / (4.0 * C * variance));
      const double beta  = (2.0 * alpha * C) / (1.0 + 2.0 * alpha * C * variance);
      
      if (alpha > 1e-10 && beta > 0.0) {
	feature_set_type::const_iterator fiter_end = features.end();
	for (feature_set_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter) {
	  const double var = covariances[fiter->first];
	  
	  weights[fiter->first]     += alpha * fiter->second * var;
	  covariances[fiter->first] -= beta * (var * var) * (fiter->second * fiter->second);
	}
      }
    }
    
    ++ samples;
  }
  
  size_t instances;
  size_t samples;
  
  double objective;
  weight_set_type weights;
  double weight_scale;
  double lambda;

  weight_set_type covariances;
};

struct OptimizerPegasos : public OptimizerBase
{
  OptimizerPegasos(const size_t& __instances,
		   const double& C) : instances(__instances), samples(0), epoch(0), weight_scale(1.0), weight_norm(0.0), lambda(C) {}
  
  void initialize()
  {
    samples = 0;
    objective = 0.0;
    
    weight_scale = 1.0;
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
  }
  void finalize()
  {
    weights *= weight_scale;
    
    weight_scale = 1.0;
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
  }

  template <typename Iterator>
  void operator()(Iterator first, Iterator last, const double loss)
  {
    const double margin = cicada::dot_product(weights, first, last, 0.0) * weight_scale;
    
    objective += loss - margin;
    
    if (loss - margin > 0.0) {
      //const double eta = 1.0 / (lambda * (epoch + 2));
      // exponential decay...
      const double eta = 0.2 * std::pow(0.85, double(epoch) / instances);
      ++ epoch;
      
      rescale(1.0 - eta * lambda);
      
      ++ samples;
      
      for (/**/; first != last; ++ first) 
	update(weights[first->first], double(first->second) * eta);
      
      // projection...
      if (weight_norm > 1.0 / lambda)
	rescale(std::sqrt(1.0 / (lambda * weight_norm)));
      
      if (weight_scale < 0.001 || 1000 < weight_scale)
	finalize();
    }
  }

  void operator()(const feature_set_type& features_reward,
		  const feature_set_type& features_penalty,
		  const double loss=1.0)
  {
    // exponential decay...
    const double eta = 0.2 * std::pow(0.85, double(epoch) / instances);
    ++ epoch;
    
    rescale(1.0 - eta * lambda);
    
    const feature_set_type features(features_reward - features_penalty);
    
    const double margin = cicada::dot_product(weights, features) * weight_scale;
    
    objective += loss - margin;
    ++ samples;
    
    feature_set_type::const_iterator fiter_end = features.end();
    for (feature_set_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
      update(weights[fiter->first], double(fiter->second) * eta);
    
    // projection...
    if (weight_norm > 1.0 / lambda)
      rescale(std::sqrt(1.0 / (lambda * weight_norm)));
    
    if (weight_scale < 0.001 || 1000 < weight_scale)
      finalize();
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
  
  double objective;
  weight_set_type weights;

  double weight_scale;
  double weight_norm;
  
  double lambda;
};


#endif
