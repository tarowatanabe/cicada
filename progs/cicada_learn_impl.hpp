//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_IMPL__HPP__
#define __CICADA_LEARN_IMPL__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/shared_ptr.hpp>

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"
#include "cicada/hypergraph.hpp"

#include "cicada/inside_outside.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/dot_product.hpp"
#include "cicada/semiring.hpp"
#include "cicada/eval.hpp"

#include "cicada_learn_online_regularize_impl.hpp"
#include "cicada_learn_online_rate_impl.hpp"

typedef boost::filesystem::path path_type;

typedef cicada::Symbol          symbol_type;
typedef cicada::Vocab           vocab_type;
typedef cicada::Sentence        sentence_type;
typedef cicada::Lattice         lattice_type;
typedef cicada::Rule            rule_type;
typedef cicada::HyperGraph      hypergraph_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef hypergraph_type::feature_set_type    feature_set_type;
typedef feature_set_type::feature_type feature_type;
typedef cicada::WeightVector<double>   weight_set_type;

struct OptimizerBase
{
  typedef OptimizerBase base_type;

  typedef cicada::semiring::Log<double> weight_type;
  typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
  typedef cicada::WeightVector<weight_type >  expectation_type;

  OptimizerBase() : samples(0), objective(0.0), weights() {}

  void initialize()
  {
    samples = 0;
    objective = 0.0;
  }

  void finalize()
  {
    
  }

  size_t samples;
  double objective;
  weight_set_type weights;
};

struct OptimizerSoftmax : public OptimizerBase
{
  OptimizerSoftmax(Regularize& __regularizer, Rate&  __rate)
    : regularizer(__regularizer.clone()),
      rate(__rate.clone()) {}
  OptimizerSoftmax(const boost::shared_ptr<Regularize>& __regularizer, const boost::shared_ptr<Rate>&  __rate)
    : regularizer(__regularizer->clone()),
      rate(__rate->clone()) {}  
  OptimizerSoftmax(const OptimizerSoftmax& x)
    : base_type(static_cast<const base_type&>(x)),
      regularizer(x.regularizer->clone()),
      rate(x.rate->clone()) {}
  OptimizerSoftmax& operator=(const OptimizerSoftmax& x)
  {
    static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
    
    regularizer.reset(x.regularizer->clone());
    rate.reset(x.rate->clone());
    
    return *this;
  }
  
  void initialize()
  {
    base_type::initialize();
    
    regularizer->initialize(weights);
  }
  
  void finalize()
  {
    base_type::finalize();
 
    regularizer->finalize(weights);
  }

  double scale() const { return regularizer->scale(); }

  gradient_type updates;
  
  void operator()(const gradient_type& correct, 
		  const gradient_type& gradient, 
		  const weight_type& Z_correct,
		  const weight_type& Z)
  {
    const double eta = rate->operator()();

    regularizer->preprocess(weights, eta);
    
    updates.clear();
    
    gradient_type::const_iterator citer_end = correct.end();
    for (gradient_type::const_iterator citer = correct.begin(); citer != citer_end; ++ citer)
      updates[citer->first] -= citer->second;
    
    gradient_type::const_iterator miter_end = gradient.end();
    for (gradient_type::const_iterator miter = gradient.begin(); miter != miter_end; ++ miter)
      updates[miter->first] += miter->second;

    gradient_type::const_iterator uiter_end = updates.end();
    for (gradient_type::const_iterator uiter = updates.begin(); uiter != uiter_end; ++ uiter) {
      const double amount = uiter->second;
      
      regularizer->update(weights, uiter->first, amount, rate->operator()(uiter->first, amount));
    }
    
    regularizer->postprocess(weights, eta);
    
    objective += double(log(Z_correct) - log(Z)) * regularizer->scale();
    
    ++ samples;
  }
  
  boost::shared_ptr<Regularize> regularizer;
  boost::shared_ptr<Rate> rate;
};

struct OptimizerMIRA : public OptimizerBase
{
  OptimizerMIRA(const double& __lambda) : lambda(__lambda) {}

  void initialize()
  {
    base_type::initialize();
  }
  
  void finalize()
  {
    base_type::finalize();
  }
  
  double scale() const { return 1.0; }

  template <typename Iterator>
  void operator()(Iterator first, Iterator last, const double loss)
  {
    const double margin   = cicada::dot_product(weights, first, last, 0.0);
    const double variance = cicada::dot_product(first, last, first, last, 0.0);
    
    objective += (loss - margin) * (loss - margin > 0.0);
    
    const double alpha = std::max(0.0, std::min(1.0 / lambda, (loss - margin) / variance));
    
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
    
    const double margin   = cicada::dot_product(weights, features);
    const double variance = cicada::dot_product(features, features);
    
    objective += (loss - margin) * (loss - margin > 0.0);
    
    const double alpha = std::max(0.0, std::min(1.0 / lambda, (loss - margin) / variance));
    
    if (alpha > 1e-10) {
      feature_set_type::const_iterator fiter_end = features.end();
      for (feature_set_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
	weights[fiter->first] += alpha * fiter->second;
    }
    ++ samples;
  }

  double lambda;
};

struct OptimizerNHERD : public OptimizerBase
{
  OptimizerNHERD(const double& __lambda) : lambda(1.0 / __lambda) {}

  void initialize()
  {
    base_type::initialize();
  }
  
  void finalize()
  {
    base_type::finalize();
  }

  double scale() const { return 1.0; }

  template <typename Iterator>
  void operator()(Iterator first, Iterator last, const double loss)
  {
    covariances.allocate(1.0);
    
    const double margin   = cicada::dot_product(weights, first, last, 0.0);
    const double variance = cicada::dot_product(first, last, covariances, first, last, 0.0); // multiply covariances...
    
    objective += (loss - margin) * (loss - margin > 0.0);
    
    const double beta = 1.0 / (variance + 1.0 / lambda);
    const double alpha = std::max(0.0, (loss - margin) * beta);
    
    if (alpha > 1e-10) {
      for (/**/; first != last; ++ first) {
	const double var = covariances[first->first];
	
	weights[first->first]     += alpha * first->second * var;
	//covariances[first->first]  = 1.0 / ((1.0 / var) + (2.0 * lambda + lambda * lambda * variance) * first->second * first->second);
	covariances[first->first]  = var / (1.0 + var * (2.0 * lambda + lambda * lambda * variance) * first->second * first->second);
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
    
    const double margin   = cicada::dot_product(weights, features);
    const double variance = cicada::dot_product(features, covariances, features); // multiply covariances...
    
    objective += (loss - margin) * (loss - margin > 0.0);
    
    const double beta = 1.0 / (variance + 1.0 / lambda);
    const double alpha = std::max(0.0, (loss - margin) * beta);
    
    if (alpha > 1e-10) {
      feature_set_type::const_iterator fiter_end = features.end();
      for (feature_set_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter) {
	const double var = covariances[fiter->first];
	
	weights[fiter->first]     += alpha * fiter->second * var;
	//covariances[fiter->first]  = 1.0 / ((1.0 / var) + (2.0 * lambda + lambda * lambda * variance) * fiter->second * fiter->second);
	covariances[fiter->first]  = var / (1.0 + var * (2.0 * lambda + lambda * lambda * variance) * fiter->second * fiter->second);
      }
    }
    
    ++ samples;
  }
  
  double lambda;
  weight_set_type covariances;
};

struct OptimizerAROW : public OptimizerBase
{
  OptimizerAROW(const double& __lambda) : lambda(__lambda) {}

  void initialize()
  {
    base_type::initialize();
  }
  
  void finalize()
  {
    base_type::finalize();
  }

  double scale() const { return 1.0; }

  template <typename Iterator>
  void operator()(Iterator first, Iterator last, const double loss)
  {
    covariances.allocate(1.0);
    
    const double margin   = cicada::dot_product(weights, first, last, 0.0);
    const double variance = cicada::dot_product(first, last, covariances, first, last, 0.0); // multiply covariances...
    
    objective += (loss - margin) * (loss - margin > 0.0);
    
    const double beta = 1.0 / (variance + lambda);
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
    
    const double margin   = cicada::dot_product(weights, features);
    const double variance = cicada::dot_product(features, covariances, features); // multiply covariances...
    
    objective += (loss - margin) * (loss - margin > 0.0);
    
    const double beta = 1.0 / (variance + lambda);
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
  
  double lambda;
  weight_set_type covariances;
};

struct OptimizerCW : public OptimizerBase
{
  OptimizerCW(const double& __lambda) : lambda(__lambda) {}

  void initialize()
  {
    base_type::initialize();
  }
  
  void finalize()
  {
    base_type::finalize();
  }

  double scale() const { return 1.0; }

  template <typename Iterator>
  void operator()(Iterator first, Iterator last, const double loss)
  {
    covariances.allocate(1.0);
    
    const double margin   = cicada::dot_product(weights, first, last, 0.0);
    const double variance = cicada::dot_product(first, last, covariances, first, last, 0.0); // multiply covariances...
    
    objective += (loss - margin) * (loss - margin > 0.0);
    
    if (loss - margin > 0.0) {
      const double theta = 1.0 + 2.0 * lambda * (margin - loss);
      const double alpha = ((- theta + std::sqrt(theta * theta - 8.0 * lambda * (margin - loss - lambda * variance))) / (4.0 * lambda * variance));
      const double beta  = (2.0 * alpha * lambda) / (1.0 + 2.0 * alpha * lambda * variance);
      
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
    
    objective += (loss - margin) * (loss - margin > 0.0);
    
    if (loss - margin > 0.0) {
      const double theta = 1.0 + 2.0 * lambda * (margin - loss);
      const double alpha = ((- theta + std::sqrt(theta * theta - 8.0 * lambda * (margin - loss - lambda * variance))) / (4.0 * lambda * variance));
      const double beta  = (2.0 * alpha * lambda) / (1.0 + 2.0 * alpha * lambda * variance);
      
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
  
  double lambda;
  weight_set_type covariances;
};

struct OptimizerHinge : public OptimizerBase
{
  OptimizerHinge(Regularize& __regularizer, Rate&  __rate)
    : regularizer(__regularizer.clone()),
      rate(__rate.clone()) {}
  OptimizerHinge(const boost::shared_ptr<Regularize>& __regularizer, const boost::shared_ptr<Rate>&  __rate)
    : regularizer(__regularizer->clone()),
      rate(__rate->clone()) {}  
  OptimizerHinge(const OptimizerHinge& x)
    : base_type(static_cast<const base_type&>(x)),
      regularizer(x.regularizer->clone()),
      rate(x.rate->clone()) {}
  OptimizerHinge& operator=(const OptimizerHinge& x)
  {
    static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
    
    regularizer.reset(x.regularizer->clone());
    rate.reset(x.rate->clone());
    
    return *this;
  }
  
  void initialize()
  {
    base_type::initialize();
    
    regularizer->initialize(weights);
  }
  void finalize()
  {
    base_type::finalize();
 
    regularizer->finalize(weights);
  }

  double scale() const { return regularizer->scale(); }

  template <typename Iterator>
  void operator()(Iterator first, Iterator last, const double loss)
  {
    const double margin = cicada::dot_product(weights, first, last, 0.0) * regularizer->scale();
    
    objective += (loss - margin) * (loss - margin > 0.0);
    
    if (loss - margin > 0.0) {
      const double eta = rate->operator()();
      
      regularizer->preprocess(weights, eta);
      
      for (/**/; first != last; ++ first) 
	regularizer->update(weights, first->first, - first->second, rate->operator()(first->first, - first->second));
      
      regularizer->postprocess(weights, eta);
      
      ++ samples;
    }
  }

  void operator()(const feature_set_type& features_reward,
		  const feature_set_type& features_penalty,
		  const double loss=1.0)
  {
    const feature_set_type features(features_reward - features_penalty);
    
    const double margin = cicada::dot_product(weights, features) * regularizer->scale();
    
    objective += (loss - margin) * (loss - margin > 0.0);
    
    const double eta = rate->operator()();
    
    regularizer->preprocess(weights, eta);
    
    feature_set_type::const_iterator fiter_end = features.end();
    for (feature_set_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
      regularizer->update(weights, fiter->first, - fiter->second, rate->operator()(fiter->first, - fiter->second));
    
    regularizer->postprocess(weights, eta);

    ++ samples;
  }

  boost::shared_ptr<Regularize> regularizer;
  boost::shared_ptr<Rate> rate;
};


#endif
