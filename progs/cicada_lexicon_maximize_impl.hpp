//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEXICON_MAXIMIZE_IMPL__HPP__
#define __CICADA_LEXICON_MAXIMIZE_IMPL__HPP__ 1

#include "utils/mathop.hpp"

struct Maximize
{
  template <typename Counts, typename Probs>
  void operator()(const Counts& counts, const Probs& probs, Probs& probs_new, const double& prior, const double& smooth)
  {
    probs_new.clear();

    double sum = 0.0;
    typename Counts::const_iterator iter_end = counts.end();
    for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter)
      sum += iter->second + prior;
    
    const double factor = 1.0 / sum;
    for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter)
      probs_new[iter->first] = (iter->second + prior) * factor;
  }
};

struct MaximizeBayes
{
  template <typename Counts, typename Probs>
  void operator()(const Counts& counts, const Probs& probs, Probs& probs_new, const double& prior, const double& smooth)
  {
    probs_new.clear();
    
    double sum = 0.0;
    typename Counts::const_iterator iter_end = counts.end();
    for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter)
      sum += iter->second + prior;
    
    const double sum_digamma = utils::mathop::digamma(sum);
    for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter)
      probs_new[iter->first] = utils::mathop::exp(utils::mathop::digamma(iter->second + prior) - sum_digamma);
  }
};

struct MaximizeL0
{
  typedef std::vector<double, std::allocator<double> > parameter_type;

  MaximizeL0(const double& __alpha, const double& __beta)
    : alpha(__alpha), beta(__beta), gamma(0.5), sigma(0.5), eta(0.1) {}

  template <typename Counts, typename Probs>
  void operator()(const Counts& counts, const Probs& probs, Probs& probs_new, const double& prior, const double& smooth)
  {
    parameter_type expected;
    parameter_type previous;
    parameter_type estimated;
    
    expected.reserve(counts.size());
    previous.reserve(counts.size());
    estimated.reserve(counts.size());
    
    if (probs.empty()) {
      // if our previous probabilities are empty (initial iteration), then
      // we perform normalization and use it as our starting point...
      
      double sum = 0.0;
      typename Counts::const_iterator iter_end = counts.end();
      for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter)
	sum += iter->second + prior;
      
      const double factor = 1.0 / sum;
      for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter) {
	expected.push_back(iter->second);
	previous.push_back((iter->second + prior) * factor);
      }
      
    } else {
      typename Counts::const_iterator iter_end = counts.end();
      for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter) {
	typename Probs::const_iterator piter = probs.find(iter->first);
	
	expected.push_back(iter->second);
	previous.push_back(piter != probs.end() ? piter->second : smooth);
      }
    }
    
    // PGD...
    gradient_descent(expected, previous, estimated);
    
    parameter_type::const_iterator piter = estimated.begin();
    typename Counts::const_iterator iter_end = counts.end();
    for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter, ++ piter)
      probs_new[iter->first] = *piter;
  }

   void gradient_descent(const parameter_type& counts, const parameter_type& point, parameter_type& point_new)
  {
    parameter_type point_curr(point);
    parameter_type point_projected(point.size());
    parameter_type point_delta(point.size());
    parameter_type gradient(point.size());
    
    for (int iter = 0; iter != 30; ++ iter) {
      const double objective_curr = compute_objective(counts, point_curr);
      
      compute_gradient(counts, point_curr, gradient);
      
      // compute new point...
      point_new.resize(counts.size());
      for (size_t i = 0; i != point_curr.size(); ++ i)
	point_new[i] = point_curr[i] - eta * gradient[i];
      
      // project into 1.0 simplex...
      project_simplex(point_new, point_projected);
      
      double armijo_bound = 0.0;
      for (size_t i = 0; i != point_curr.size(); ++ i)
	armijo_bound += sigma * gamma * gradient[i] * (point_projected[i] - point_curr[i]);
      
      bool updated = false;
      
      double gamma_curr = gamma;
      double gamma_min = 0.0;
      double armijo_bound_curr = armijo_bound;
      double objective_min = objective_curr;
      
      for (int steps = 0; steps != 20; ++ steps) {
	for (size_t i = 0; i != point_curr.size(); ++ i)
	  point_delta[i] = (1.0 - gamma_curr) * point_curr[i] + gamma_curr * point_projected[i];
	
	const double objective_delta = compute_objective(counts, point_delta);
	
	if (objective_delta < objective_min) {
	  objective_min = objective_delta;
	  gamma_min = gamma_curr;
	  updated = true;
	}
	
	if (objective_delta <= objective_curr + armijo_bound_curr)
	  break;
	
	gamma_curr *= gamma;
	armijo_bound_curr *= gamma;
      }
      
      // finish gradient descent, since we made no update!
      if (! updated) break;
      
      // update current point...
      for (size_t i = 0; i != point_curr.size(); ++ i)
	point_curr[i] = (1.0 - gamma_min) * point_curr[i] + gamma_min * point_projected[i];
    }
    
    // the final current point == result..
    point_new = point_curr;
  }

  double compute_objective(const parameter_type& counts, const parameter_type& point) const
  {
    parameter_type::const_iterator citer_end = counts.end();
    parameter_type::const_iterator citer     = counts.begin();
    
    parameter_type::const_iterator piter_end = point.end();
    parameter_type::const_iterator piter     = point.begin();

    double objective = 0.0;
    
    for (/**/; citer != citer_end; ++ citer, ++ piter) {
      if (*piter != 0.0)
	objective -= *citer * std::log(*piter) + alpha * std::exp(- *piter / beta);
      else
	objective -= alpha * std::exp(- *piter / beta);
    }

    return objective;
  }
  
  
  void compute_gradient(const parameter_type& counts, const parameter_type& point, parameter_type& gradient) const
  {
    gradient.resize(point.size());

    parameter_type::const_iterator citer_end = counts.end();
    parameter_type::const_iterator citer     = counts.begin();
    
    parameter_type::const_iterator piter_end = point.end();
    parameter_type::const_iterator piter     = point.begin();
    
    parameter_type::iterator giter_end = gradient.end();
    parameter_type::iterator giter = gradient.begin();
    
    for (/**/; citer != citer_end; ++ citer, ++ piter, ++ giter) {
      if (*piter != 0.0)
	*giter = - (*citer / *piter) + alpha / beta * std::exp(- *piter / beta);
      else
	*giter = alpha / beta * std::exp(- *piter / beta);
    }
  }
  
  //
  // projection onto a simplex
  //
  // @inproceedings{Duchi:2008:EPL:1390156.1390191,
  //  author = {Duchi, John and Shalev-Shwartz, Shai and Singer, Yoram and Chandra, Tushar},
  //  title = {Efficient projections onto the l1-ball for learning in high dimensions},
  //  booktitle = {Proceedings of the 25th international conference on Machine learning},
  //  series = {ICML '08},
  //  year = {2008},
  //  isbn = {978-1-60558-205-4},
  //  location = {Helsinki, Finland},
  //  pages = {272--279},
  //  numpages = {8},
  //  url = {http://doi.acm.org/10.1145/1390156.1390191},
  //  doi = {10.1145/1390156.1390191},
  //  acmid = {1390191},
  //  publisher = {ACM},
  //  address = {New York, NY, USA},
  // } 
  
  struct projection
  {
    projection(const double& __theta) : theta(__theta) {}

    double operator()(const double& x) const
    {
      return std::max(x - theta, 0.0);
    }
    
    double theta;
  };
  
  void project_simplex(const parameter_type& v, parameter_type& projected)
  {
    if (v.empty()) {
      projected.clear();
      return;
    }
    
    mu.clear();
    mu.insert(mu.end(), v.begin(), v.end());
    std::sort(mu.begin(), mu.end(), std::greater<double>());
    
    const double z = 1.0;
    
    size_t j = 1;
    double sum = 0.0;
    
    size_t rho = 0;
    double rho_sum = 0.0;
    
    parameter_type::const_iterator miter_end = mu.end();
    for (parameter_type::const_iterator miter = mu.begin(); miter != miter_end; ++ miter, ++ j) {
      sum += *miter;
      
      if (*miter - (1.0 / j) * (sum - z) > 0.0) {
	rho = j;
	rho_sum = sum;
      }
    }
    
    const double theta = (1.0 / rho) * (rho_sum - z);
    
    projected.resize(v.size());
    
    std::transform(v.begin(), v.end(), projected.begin(), projection(theta));
  }
  
  double alpha;
  double beta;
  
  double gamma;
  double sigma;
  double eta;
  
  parameter_type mu;
};

#endif
