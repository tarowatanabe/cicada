//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_ONLINE_REGULARIZE_IMPL__HPP__
#define __CICADA_LEARN_ONLINE_REGULARIZE_IMPL__HPP__ 1

#include <numeric>
#include <vector>
#include <deque>
#include <cmath>

#include <cicada/hypergraph.hpp>
#include <cicada/feature_vector.hpp>
#include <cicada/weight_vector.hpp>

#include <utils/mathop.hpp>

struct RegularizeL1
{
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

struct RegularizeOSCAR
{
  typedef cicada::HyperGraph hypergraph_type;
  typedef hypergraph_type::feature_set_type    feature_set_type;
  
  typedef feature_set_type::feature_type feature_type;

  typedef cicada::WeightVector<double> weight_set_type;
  
  // TODO: fix me!
  // I hate this constant! 
  
  RegularizeOSCAR(const double lambda1, const double lambda2=oscar) : lambda1_(lambda1), lambda2_(lambda2) {}
  
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
  typedef std::vector<id_type, std::allocator<id_type> > index_set_type;

  typedef std::pair<id_type, id_type> item_type;
  typedef std::vector<item_type, std::allocator<item_type> > group_type;
  typedef std::pair<group_type, double> group_value_type;

  typedef std::deque<group_value_type, std::allocator<group_value_type> > stack_type;

  index_set_type p;
  stack_type     G;
  group_type     g;

  struct greater_weights
  {
    greater_weights(const weight_set_type& weights) : weights_(weights) {}

    bool operator()(const id_type& x, const id_type& y) const
    {
      return std::fabs(weights_[x]) > std::fabs(weights_[y]);
    }
    
    const weight_set_type& weights_;
  };

  weight_set_type weights_new;
    
public:
  void postprocess(weight_set_type& weights, const double& rate)
  {
    p.clear();
    for (id_type id = 0; id != weights.size(); ++ id)
      if (weights[id] != 0.0)
	p.push_back(id);
    
    // nothing to do!
    if (p.empty()) return;
    
    // sort...
    std::sort(p.begin(), p.end(), greater_weights(weights));

    // initialize stack...
    g.clear();
    g.push_back(std::make_pair(1, p.front()));

    G.clear();
    G.push_back(group_value_type(g, objective(weights, g, rate)));
    
    // iterate p...
    for (id_type i = 2; i != p.size() + 1; ++ i) {
      g.clear();
      g.push_back(std::make_pair(i, p[i]));

      double value = objective(weights, g, rate);
      
      while (! G.empty() && value >= G.back().second) {
	g.insert(g.end(), G.back().first.begin(), G.back().first.end());
	value = objective(weights, g, rate);
	
	G.pop_back();
      }
      
      G.push_back(group_value_type(g, value));
    }
    
    weights_new.clear();
    
    stack_type::const_iterator giter_end = G.end();
    for (stack_type::const_iterator giter = G.begin(); giter != giter_end; ++ giter) {
      const double& value = giter->second;
      
      group_type::const_iterator iter_end = giter->first.end();
      for (group_type::const_iterator iter = giter->first.begin(); iter != iter_end; ++ iter)
	weights_new[iter->second] = utils::mathop::sgn(weights[iter->second]) * value;
    }
    
    weights.swap(weights_new);
  }
  
  double objective(const weight_set_type& weights, const group_type& group, const double& rate) const
  {
    double v = 0.0;
    
    group_type::const_iterator giter_end = group.end();
    for (group_type::const_iterator giter = group.begin(); giter != giter_end; ++ giter)
      v += std::fabs(weights[giter->second]) - 2.0 * rate * (lambda1_ + lambda2_ * (weights.size() - giter->first));
    
    return 0.5 * v / group.size();
  }

private:
  
  double lambda1_;
  double lambda2_;
};

#endif
