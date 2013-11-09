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

class Regularize
{
public:
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef cicada::HyperGraph hypergraph_type;
  typedef hypergraph_type::feature_set_type    feature_set_type;
  
  typedef feature_set_type::feature_type feature_type;

  typedef cicada::WeightVector<double> weight_set_type;

public:
  Regularize() {}
  virtual ~Regularize() {}
  
  virtual
  Regularize* clone() const = 0;

  virtual 
  double scale() const = 0;
  
  virtual
  void initialize(weight_set_type& weights) = 0;
  
  virtual
  void finalize(weight_set_type& weights) = 0;

  virtual
  void preprocess(weight_set_type& weights, const double& eta) = 0;

  virtual
  void update(weight_set_type& weights, const feature_type& feature, const double& amount, const double& eta) = 0;
  
  virtual
  void postprocess(weight_set_type& weights, const double& eta) = 0;
};

class RegularizeNone : public Regularize
{
public:
  Regularize* clone() const { return new RegularizeNone(*this); }

  double scale() const { return 1.0; }

  void initialize(weight_set_type& weights)
  {
    
  }
  
  void finalize(weight_set_type& weights)
  {
    
  }

  void preprocess(weight_set_type& weights, const double& eta)
  {

  }

  void update(weight_set_type& weights, const feature_type& feature, const double& amount, const double& eta)
  {
    weights[feature] -= amount * eta;
  }

  void postprocess(weight_set_type& weights, const double& eta)
  {
    
  }
};

// L1 regularization from
//
// @InProceedings{tsuruoka-tsujii-ananiadou:2009:ACLIJCNLP,
//   author    = {Tsuruoka, Yoshimasa  and  Tsujii, Jun'ichi  and  Ananiadou, Sophia},
//   title     = {Stochastic Gradient Descent Training for L1-regularized Log-linear Models with Cumulative Penalty},
//   booktitle = {Proceedings of the Joint Conference of the 47th Annual Meeting of the ACL and the 4th International Joint Conference on Natural Language Processing of the AFNLP},
//   month     = {August},
//   year      = {2009},
//   address   = {Suntec, Singapore},
//   publisher = {Association for Computational Linguistics},
//   pages     = {477--485},
//   url       = {http://www.aclweb.org/anthology/P/P09/P09-1054}
// }

class RegularizeL1 : public Regularize
{
public:
  typedef  weight_set_type penalty_set_type;

  RegularizeL1(const double lambda)
    : lambda_(lambda), penalties_(), penalty_(0.0) {}

  Regularize* clone() const { return new RegularizeL1(*this); }

  double scale() const { return 1.0; }
  
  void initialize(weight_set_type& weights)
  {
    
  }
  
  void finalize(weight_set_type& weights)
  {
    
  }
  
  void preprocess(weight_set_type& weights, const double& eta)
  {
    //penalty_ += eta * lambda_;
  }

  void update(weight_set_type& weights, const feature_type& feature, const double& amount, const double& eta)
  {
    double& x = weights[feature];
    //double& penalty = penalties_[feature];
    
    // we use FOBOS, which is far simpler...
    
    const double f = x - eta * amount;
    
    x = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - eta * lambda_);

#if 0
    // update...
    x -= amount * eta;
    
    // apply penalties...
    const double x_half = x;
    
    if (x > 0.0)
      x = std::max(0.0, x - penalty - penalty_);
    else if (x < 0.0)
      x = std::min(0.0, x - penalty + penalty_);
    
    penalty += x - x_half;
#endif
  }

  void postprocess(weight_set_type& weights, const double& eta)
  {
    
  }

private:  
  double lambda_;
  
  penalty_set_type penalties_;
  double penalty_;
};


// L1 and L2 regularization, based on the code from L1 regularizer
class RegularizeL1L2 : public Regularize
{
public:
  typedef  weight_set_type penalty_set_type;

  RegularizeL1L2(const double lambda1, const double lambda2)
    : lambda1_(lambda1), lambda2_(lambda2), scale_(1.0), penalties_(), penalty_(0.0) {}

  Regularize* clone() const { return new RegularizeL1L2(*this); }
  
  double scale() const { return scale_; }
  
  void initialize(weight_set_type& weights)
  {
    scale_ = 1.0;
  }
  
  void finalize(weight_set_type& weights)
  {
    if (scale_ != 1.0)
      weights *= scale_;
    scale_ = 1.0;
  }
  
  void preprocess(weight_set_type& weights, const double& eta)
  {
    // rescaling weights!
    scale_ *= 1.0 - eta * lambda2_;
    if (scale_ == 0.0) {
      scale_ = 1.0;
      std::fill(weights.begin(), weights.end(), 0.0);
    }
    
    //penalty_ += eta * lambda1_;
  }

  void update(weight_set_type& weights, const feature_type& feature, const double& amount, const double& eta)
  {
    double& x = weights[feature];
    
    const double f = x * scale_ - eta * amount;
    
    x = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - eta * lambda1_) / scale_;
    
    //weights[feature] -= amount * eta / scale_;
  }
  
  void postprocess(weight_set_type& weights, const double& eta)
  {
#if 0
    typedef feature_type::id_type id_type;
    
    const double factor = 1.0 / scale_;
    
    for (id_type id = 0; id != weights.size(); ++ id)
      if (weights[id] != 0.0) {
	double& x = weights[id];
	double& penalty = penalties_[id];
	
	// apply penalties...
	const double x_half = x;
	
	if (x > 0.0)
	  x = std::max(0.0, x - (penalty + penalty_) * factor);
	else if (x < 0.0)
	  x = std::min(0.0, x - (penalty - penalty_) * factor);
	
	penalty += (x - x_half) * scale_;
      }
#endif

    if (scale_ < 0.001 || 1000 < scale_)
      finalize(weights);
  }

private:  
  double lambda1_;
  double lambda2_;

  double scale_;
  
  penalty_set_type penalties_;
  double penalty_;
};

// L2 regularization from:
//
// @inproceedings{Shalev-Shwartz:2007:PPE:1273496.1273598,
//  author = {Shalev-Shwartz, Shai and Singer, Yoram and Srebro, Nathan},
//  title = {Pegasos: Primal Estimated sub-GrAdient SOlver for SVM},
//  booktitle = {Proceedings of the 24th international conference on Machine learning},
//  series = {ICML '07},
//  year = {2007},
//  isbn = {978-1-59593-793-3},
//  location = {Corvalis, Oregon},
//  pages = {807--814},
//  numpages = {8},
//  url = {http://doi.acm.org/10.1145/1273496.1273598},
//  doi = {10.1145/1273496.1273598},
//  acmid = {1273598},
//  publisher = {ACM},
//  address = {New York, NY, USA},
// } 

class RegularizeL2 : public Regularize
{
public:
  RegularizeL2(const double lambda) : lambda_(lambda), scale_(1.0), norm_(0.0) {}  

  Regularize* clone() const { return new RegularizeL2(*this); }
  
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

  void preprocess(weight_set_type& weights, const double& eta)
  {
    rescale(weights, 1.0 - eta * lambda_);
  }
  
  void update(weight_set_type& weights, const feature_type& feature, const double& __amount, const double& eta)
  {
    double& x = weights[feature];

    const double amount = - __amount * eta;
    
    norm_ += 2.0 * x * scale_ * amount  + amount * amount;
    
    x += amount / scale_;
  }
  
  void postprocess(weight_set_type& weights, const double& eta)
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
//
//
// @ARTICLE{6238378, 
// author={Zhong, L.W. and Kwok, J.T.}, 
// journal={Neural Networks and Learning Systems, IEEE Transactions on}, 
// title={Efficient Sparse Modeling With Automatic Feature Grouping}, 
// year={2012}, 
// volume={23}, 
// number={9}, 
// pages={1436-1447}, 
// keywords={computational complexity;computational geometry;data handling;gradient methods;learning (artificial intelligence);pattern clustering;regression analysis;OSCAR;accelerated gradient method;automatic feature grouping;data interpretation;data understanding;empirical time complexity reduction;estimation variance reduction;feature coefficients;feature selection stability improvement;high-dimensional data;iterative group merging algorithm;l1 -regularizer;learning process;octagonal shrinkage-and-clustering algorithm-for-regression;optimization procedure;pairwise l∞-regularizer;sparse modeling;sparse-modeling approach;Acceleration;Algorithm design and analysis;Complexity theory;Computational modeling;Convergence;Gradient methods;Accelerated gradient descent;feature grouping;sparse modeling;structured sparsity}, 
// doi={10.1109/TNNLS.2012.2200262}, 
// ISSN={2162-237X},}

class RegularizeOSCAR : public Regularize
{
public:
  
  RegularizeOSCAR(const double lambda1, const double lambda2) : lambda1_(lambda1), lambda2_(lambda2) {}

  Regularize* clone() const { return new RegularizeOSCAR(*this); }
  
  double scale() const { return 1.0; }

  void initialize(weight_set_type& weights)
  {
    
  }
  
  void finalize(weight_set_type& weights)
  {
    
  }

  void preprocess(weight_set_type& weights, const double& eta)
  {
    
  }
  
  void update(weight_set_type& weights, const feature_type& feature, const double& amount, const double& eta)
  {
    weights[feature] -= amount * eta;
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
    group_type(const size_type& first, const double& score) : first_(first), last_(first + 1), score_(score) {}
    
    double score() const { return score_ / (last_ - first_); }
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
    G.push_back(group_type(0, std::fabs(p.front().second) - eta * (lambda1_ + lambda2_ * (num_features - 1))));
    
    // iterate p and perform grouping...
    for (id_type i = 1; i != p.size(); ++ i) {
      group_type g(i, std::fabs(p[i].second) - eta * (lambda1_ + lambda2_ * (num_features - i - 1)));
      
      while (! G.empty() && g.score() >= G.back().score()) {
	// merge group
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
	const double score = giter->score();
	
	for (size_type i = giter->first_; i != giter->last_; ++ i)
	  weights[p[i].first] = utils::mathop::sgn(p[i].second) * score;
      }
  }

private:
  
  double lambda1_;
  double lambda2_;
};

// Regularized Dual Averaging variant of regularization from:
//
// @article{Xiao:2010:DAM:1756006.1953017,
//  author = {Xiao, Lin},
//  title = {Dual Averaging Methods for Regularized Stochastic Learning and Online Optimization},
//  journal = {J. Mach. Learn. Res.},
//  issue_date = {3/1/2010},
//  volume = {11},
//  month = dec,
//  year = {2010},
//  issn = {1532-4435},
//  pages = {2543--2596},
//  numpages = {54},
//  url = {http://dl.acm.org/citation.cfm?id=1756006.1953017},
//  acmid = {1953017},
//  publisher = {JMLR.org},
// } 

class RegularizeRDAL1 : public Regularize
{
public:

  RegularizeRDAL1(const double lambda)
    : lambda_(lambda), averaged_(), scale_(1.0), epoch_(0) {}

  Regularize* clone() const { return new RegularizeRDAL1(*this); }
  
  double scale() const { return 1.0; }
  
  void initialize(weight_set_type& weights)
  {
    
  }
  
  void finalize(weight_set_type& weights)
  {
    
  }
  
  void preprocess(weight_set_type& weights, const double& eta)
  {
    if (epoch_)
      scale_ *= double(epoch_) / (epoch_ + 1);
  }

  void update(weight_set_type& weights, const feature_type& feature, const double& amount, const double& eta)
  {
    if (epoch_)
      averaged_[feature] += amount / (scale_ * (epoch_ + 1));
    else
      averaged_[feature] = amount;
  }
  
  void postprocess(weight_set_type& weights, const double& eta)
  {
    typedef feature_type::id_type id_type;
    
    weights.clear();

    const double lambda = lambda_;
    
    for (id_type id = 0; id != averaged_.size(); ++ id) {
      const double value = averaged_[id] * scale_;

      if (value < - lambda)
	weights[feature_type(id)] = - (value + lambda);
      else if (value > lambda)
	weights[feature_type(id)] = - (value - lambda);
    }
    
    ++ epoch_;
  }
  
private:  
  double lambda_;

  weight_set_type averaged_;
  double scale_;
  size_type epoch_;
};

class RegularizeRDAL1L2 : public Regularize
{
public:
  RegularizeRDAL1L2(const double lambda1, const double lambda2)
    : weight_scale_(1.0), lambda1_(lambda1), lambda2_(lambda2), averaged_(), scale_(1.0), epoch_(0) {}

  Regularize* clone() const { return new RegularizeRDAL1L2(*this); }
  
  double scale() const { return weight_scale_; }
  
  void initialize(weight_set_type& weights)
  {
    weight_scale_ = 1.0;
  }
  
  void finalize(weight_set_type& weights)
  {
    if (weight_scale_ != 1.0)
      weights *= weight_scale_;
    weight_scale_ = 1.0;
  }
  
  void preprocess(weight_set_type& weights, const double& eta)
  {
    if (epoch_)
      scale_ *= double(epoch_) / (epoch_ + 1);
  }

  void update(weight_set_type& weights, const feature_type& feature, const double& amount, const double& eta)
  {
    if (epoch_)
      averaged_[feature] += amount / (scale_ * (epoch_ + 1));
    else
      averaged_[feature] = amount;
  }
  
  void postprocess(weight_set_type& weights, const double& eta)
  {
    typedef feature_type::id_type id_type;
    
    weights.clear();

    const double lambda = lambda1_;
    
    double norm = 0.0;

    for (id_type id = 0; id != averaged_.size(); ++ id) {
      const double value = averaged_[id] * scale_;

      if (value < - lambda) {
	weights[feature_type(id)] = - (value + lambda);
	
	norm += (value + lambda) * (value + lambda);
      } else if (value > lambda) {
	weights[feature_type(id)] = - (value - lambda);

	norm += (value - lambda) * (value - lambda);
      }
    }

    weight_scale_ = 1.0;
    if (norm > 1.0 / lambda2_) {
      weight_scale_ = std::sqrt((1.0 / lambda2_) * (1.0 / norm));
      
      if (weight_scale_ < 0.001 || 1000 < weight_scale_)
	finalize(weights);
    }
    
    ++ epoch_;
  }
  
private:  
  double weight_scale_;

  double lambda1_;
  double lambda2_;

  weight_set_type averaged_;
  double scale_;
  size_type epoch_;
};

class RegularizeRDAL2 : public Regularize
{
public:
  RegularizeRDAL2(const double lambda)
    : weight_scale_(1.0), lambda_(lambda), averaged_(), scale_(1.0), epoch_(0) {}

  Regularize* clone() const { return new RegularizeRDAL2(*this); }
  
  double scale() const { return weight_scale_; }
  
  void initialize(weight_set_type& weights)
  {
    weight_scale_ = 1.0;
  }
  
  void finalize(weight_set_type& weights)
  {
    if (weight_scale_ != 1.0)
      weights *= weight_scale_;
    weight_scale_ = 1.0;
  }
  
  void preprocess(weight_set_type& weights, const double& eta)
  {
    if (epoch_)
      scale_ *= double(epoch_) / (epoch_ + 1);
  }

  void update(weight_set_type& weights, const feature_type& feature, const double& amount, const double& eta)
  {
    if (epoch_)
      averaged_[feature] += amount / (scale_ * (epoch_ + 1));
    else
      averaged_[feature] = amount;
  }
  
  void postprocess(weight_set_type& weights, const double& eta)
  {
    typedef feature_type::id_type id_type;
    
    weights.clear();
    
    double norm = 0.0;
    
    for (id_type id = 0; id != averaged_.size(); ++ id) {
      const double value = averaged_[id] * scale_;
      
      weights[feature_type(id)] = - value;
      norm += value * value;
    }
    
    weight_scale_ = 1.0;
    if (norm > 1.0 / lambda_) {
      weight_scale_ = std::sqrt((1.0 / lambda_) * (1.0 / norm));
      
      if (weight_scale_ < 0.001 || 1000 < weight_scale_)
	finalize(weights);
    }
    
    ++ epoch_;
  }
  
private:  
  double weight_scale_;
  
  double lambda_;

  weight_set_type averaged_;
  double scale_;
  size_type epoch_;
};

class RegularizeRDAOSCAR : public Regularize
{
public:
  RegularizeRDAOSCAR(const double lambda1, const double lambda2)
    : lambda1_(lambda1), lambda2_(lambda2), averaged_(), scale_(1.0), epoch_(0) {}
  
  Regularize* clone() const { return new RegularizeRDAOSCAR(*this); }
  
  double scale() const { return 1.0; }
  
  void initialize(weight_set_type& weights)
  {
    
  }
  
  void finalize(weight_set_type& weights)
  {
    
  }
  
  void preprocess(weight_set_type& weights, const double& eta)
  {
    if (epoch_)
      scale_ *= double(epoch_) / (epoch_ + 1);
  }

  void update(weight_set_type& weights, const feature_type& feature, const double& amount, const double& eta)
  {
    if (epoch_)
      averaged_[feature] += amount / (scale_ * (epoch_ + 1));
    else
      averaged_[feature] = amount;
  }
  
private:
  typedef std::pair<feature_type, double> index_type;
  typedef std::vector<index_type, std::allocator<index_type> > index_set_type;

  struct group_type
  {
    size_type first_;
    size_type last_;
    double    score_;

    group_type() : first_(), last_(), score_() {}
    group_type(const size_type& first, const double& score) : first_(first), last_(first + 1), score_(score) {}
    
    double score() const { return score_ / (last_ - first_); }
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
    typedef feature_type::id_type id_type;

    weights.clear();
    
    const size_type num_features = feature_type::allocated();
    
    p.clear();
    for (id_type id = 0; id != averaged_.size(); ++ id)
      if (averaged_[id] != 0.0)
	p.push_back(std::make_pair(feature_type(id), - averaged_[id] * scale_));
    
    // nothing to do!
    if (p.empty()) return;
    
    // sort...
    std::sort(p.begin(), p.end(), greater_weights());
    
    // initialize stack...
    G.clear();
    G.push_back(group_type(0, std::fabs(p.front().second) - (lambda1_ + lambda2_ * (num_features - 1))));
    
    // iterate p and perform grouping...
    for (id_type i = 1; i != p.size(); ++ i) {
      group_type g(i, std::fabs(p[i].second) - (lambda1_ + lambda2_ * (num_features - i - 1)));
      
      while (! G.empty() && g.score() >= G.back().score()) {
	// merge group
	g.first_ = G.back().first_;
	g.score_ += G.back().score_;
	
	G.pop_back();
      }
      
      G.push_back(g);
    }
    
    stack_type::const_iterator giter_end = G.end();
    for (stack_type::const_iterator giter = G.begin(); giter != giter_end; ++ giter) 
      if (giter->score_ > 0.0) {
	const double score = giter->score();
	
	for (size_type i = giter->first_; i != giter->last_; ++ i)
	  weights[p[i].first] = utils::mathop::sgn(p[i].second) * score;
      }
    
    ++ epoch_;
  }
  
private:  
  double lambda1_;
  double lambda2_;

  weight_set_type averaged_;
  double scale_;
  size_type epoch_;
};

#endif
