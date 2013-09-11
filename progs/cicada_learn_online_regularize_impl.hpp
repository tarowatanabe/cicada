//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_ONLINE_REGULARIZE_IMPL__HPP__
#define __CICADA_LEARN_ONLINE_REGULARIZE_IMPL__HPP__ 1

#include <numeric>
#include <vector>
#include <cmath>

#include <cicada/hypergraph.hpp>
#include <cicada/feature_vector.hpp>
#include <cicada/weight_vector.hpp>

#include <utils/mathop.hpp>

struct RegularizeL1
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef cicada::HyperGraph hypergraph_type;
  typedef hypergraph_type::feature_set_type    feature_set_type;
  
  typedef feature_set_type::feature_type feature_type;

  typedef cicada::WeightVector<double> penalty_set_type;
  typedef cicada::WeightVector<double> weight_set_type;

  RegularizeL1(const double lambda) : lambda_(lambda), penalties_(), penalty_(0.0) {}

  double scale() const { return 1.0; }
  
  void initialize(weight_set_type& weights)
  {
    
  }
  
  void finalize(weight_set_type& weights)
  {
    
  }
  
  void preprocess(weight_set_type& weights, const double& rate)
  {
    penalty_ += rate * lambda_;
  }

  void update(weight_set_type& weights, const feature_type& feature, const double amount)
  {
    double& x = weights[feature];
    double& penalty = penalties_[feature];
    
    // update...
    x += amount;
    
    // apply penalties...
    const double x_half = x;
    
    if (x > 0.0)
      x = std::max(0.0, x - penalty - penalty_);
    else if (x < 0.0)
      x = std::min(0.0, x - penalty + penalty_);
    
    penalty += x - x_half;
  }

  void postprocess(weight_set_type& weights, const double& rate)
  {
    
  }

private:  
  double lambda_;
  
  penalty_set_type penalties_;
  double penalty_;
};

struct RegularizeL2
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef cicada::HyperGraph hypergraph_type;
  typedef hypergraph_type::feature_set_type    feature_set_type;
  
  typedef feature_set_type::feature_type feature_type;

  typedef cicada::WeightVector<double> weight_set_type;

  RegularizeL2(const double lambda) : lambda_(lambda), scale_(1.0), norm_(0.0) {}  
  
  double scale() const { return scale_; }
  
  void initialize(weight_set_type& weights)
  {
    scale_ = 1.0;
    norm_ = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
  }
  
  void finalize(weight_set_type& weights)
  {
    weights *= scale_;
    
    scale_ = 1.0;
    norm_ = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
  }

  void preprocess(weight_set_type& weights, const double& rate)
  {
    rescale(weights, 1.0 - rate * lambda_);
  }
  
  void update(weight_set_type& weights, const feature_type& feature, const double amount)
  {
    double& x = weights[feature];
    
    norm_ += 2.0 * x * scale_ * amount  + amount * amount;
    
    x += amount / scale_;
  }
  
  void postprocess(weight_set_type& weights, const double& rate)
  {
    if (norm_ > 1.0 / lambda_)
      rescale(weights, std::sqrt((1.0 / lambda_) * (1.0 / norm_)));
    
    if (scale_ < 0.001 || 1000 < scale_)
      finalize(weights);
  }
  
private:
  void rescale(weight_set_type& weights, const double scaling)
  {
    norm_ *= scaling * scaling;
    
    if (scaling != 0.0)
      scale_ *= scaling;
    else {
      scale_ = 1.0;
      std::fill(weights.begin(), weights.end(), 0.0);
    }
  }

private:  
  double lambda_;

  double scale_;
  double norm_;
};

//
// an online version of OSCAR...
//
// Note that the original OSCAR has factor of two, which can be ignored in the computation.
struct RegularizeOSCAR
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef cicada::HyperGraph hypergraph_type;
  typedef hypergraph_type::feature_set_type    feature_set_type;
  
  typedef feature_set_type::feature_type feature_type;

  typedef cicada::WeightVector<double> weight_set_type;
  
  // TODO: fix me!
  // I hate this constant! 
  
  RegularizeOSCAR(const double lambda1, const double lambda2) : lambda1_(lambda1), lambda2_(lambda2) {}
  
  double scale() const { return 1.0; }

  void initialize(weight_set_type& weights)
  {
    
  }
  
  void finalize(weight_set_type& weights)
  {
    
  }

  void preprocess(weight_set_type& weights, const double rate)
  {
    
  }
  
  void update(weight_set_type& weights, const feature_type& feature, const double amount)
  {
    weights[feature] += amount;
  }

private:
  typedef feature_type::id_type id_type;
  typedef std::pair<feature_type, double> index_type;
  typedef std::vector<index_type, std::allocator<index_type> > index_set_type;

  struct group_type
  {
    size_type first_;
    size_type last_;
    double    score_;

    group_type() : first_(), last_(), score_() {}
    group_type(const size_type& first, const size_type& last, const double& score) : first_(first), last_(last), score_(score) {}
    
    size_type size() const { return last_ - first_; }
  };
  
  typedef std::vector<group_type, std::allocator<group_type> > stack_type;
  
  index_set_type p;
  stack_type     G;

  struct greater_weights
  {
    bool operator()(const index_type& x, const index_type& y) const
    {
      return std::fabs(x.second) > std::fabs(y.second);
    }
  };

public:
  void postprocess(weight_set_type& weights, const double& eta)
  {
    const size_type num_features = feature_type::allocated();
    
    p.clear();
    for (id_type id = 0; id != weights.size(); ++ id)
      if (weights[id] != 0.0)
	p.push_back(std::make_pair(feature_type(id), weights[id]));
    
    // nothing to do!
    if (p.empty()) return;
    
    // sort...
    std::sort(p.begin(), p.end(), greater_weights());
    
    // initialize stack...
    G.clear();
    G.push_back(group_type(0, 1, std::fabs(p.front().second) - eta * (lambda1_ + lambda2_ * (num_features - 1))));
    
    // iterate p and perform grouping...
    for (id_type i = 1; i != p.size(); ++ i) {
      group_type g(i, i + 1, std::fabs(p[i].second) - eta * (lambda1_ + lambda2_ * (num_features - i - 1)));
      
      while (! G.empty() && (g.score_ / g.size()) >= (G.back().score_ / G.back().size())) {
	g.first_ = G.back().first_;
	g.score_ += G.back().score_;
	
	G.pop_back();
      }
      
      G.push_back(g);
    }
    
    // compute new parameters
    weights.clear();
    
    stack_type::const_iterator giter_end = G.end();
    for (stack_type::const_iterator giter = G.begin(); giter != giter_end; ++ giter) 
      if (giter->score_ > 0.0) {
	const double score = giter->score_ / giter->size();
	
	for (size_type i = giter->first_; i != giter->last_; ++ i)
	  weights[p[i].first] = utils::mathop::sgn(p[i].second) * score;
      }
  }

private:
  
  double lambda1_;
  double lambda2_;
};

#endif
