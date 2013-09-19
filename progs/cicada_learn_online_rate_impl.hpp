//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_ONLINE_RATE_IMPL__HPP__
#define __CICADA_LEARN_ONLINE_RATE_IMPL__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/weight_vector.hpp>

#include <utils/mathop.hpp>

class Rate
{
public:
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
 
  typedef cicada::HyperGraph hypergraph_type;
  typedef hypergraph_type::feature_set_type    feature_set_type;
  
  typedef feature_set_type::feature_type feature_type;
  
public:
  Rate() {}
  virtual ~Rate() {}

  virtual Rate* clone() const = 0;
  
  virtual double operator()() const = 0;

  virtual double operator()(const feature_type& feat, const double& grad) const = 0;
};

class RateNone : public Rate
{
public:
  
  Rate* clone() const { return new RateNone(*this); }

  double operator()() const
  {
    return 1.0;
  }

  double operator()(const feature_type& feat, const double& grad) const
  {
    return 1.0;
  }
};

class RateSimple : public Rate
{
  double eta0_;
  
  double eta_;
  size_type epoch_;

public:  
  RateSimple(const double& eta0) : eta0_(eta0), epoch_(0) {}
  
  Rate* clone() const { return new RateSimple(*this); }

  double operator()() const
  {
    const_cast<double&>(eta_) = eta0_ / (epoch_ + 1);
    ++ const_cast<size_type&>(epoch_);
    return eta_;
  }
  
  double operator()(const feature_type& feat, const double& grad) const
  {
    return eta_;
  }
};

class RateSqrt : public Rate
{
  double eta0_;
  
  double eta_;
  size_type epoch_;

public:  
  RateSqrt(const double& eta0) : eta0_(eta0), epoch_(0) {}

  Rate* clone() const { return new RateSqrt(*this); }
  
  double operator()() const
  {
    const_cast<double&>(eta_) = eta0_ / std::sqrt(double(epoch_ + 1));
    ++ const_cast<size_type&>(epoch_);
    return eta_;
  }
  
  double operator()(const feature_type& feat, const double& grad) const
  {
    return eta_;
  }  
};

class RateExponential : public Rate
{
  double alpha0_;
  double eta0_;
  size_type samples0_;
  
  double eta_;
  size_type epoch_;

public:  
  RateExponential(const double& alpha0, const double& eta0, const size_type& samples0) : alpha0_(alpha0), eta0_(eta0), samples0_(samples0), epoch_(0) {}

  Rate* clone() const { return new RateExponential(*this); }
  
  double operator()() const
  {
    const_cast<double&>(eta_) = eta0_ * std::pow(alpha0_, double(epoch_) / samples0_);
    ++ const_cast<size_type&>(epoch_);
    return eta_;
  }
  
  double operator()(const feature_type& feat, const double& grad) const
  {
    return eta_;
  }
};

class RateAdaGrad : public Rate
{
  typedef cicada::WeightVector<double> weight_set_type;

  double eta0_;

  weight_set_type grads2_;
  double eta_;
  size_type epoch_;

public:  
  RateAdaGrad(const double& eta0) : eta0_(eta0), epoch_(0) {}

  Rate* clone() const { return new RateAdaGrad(*this); }
  
  double operator()() const
  {
    const_cast<double&>(eta_) = eta0_ / (epoch_ + 1);
    ++ const_cast<size_type&>(epoch_);
    return eta_;
  }
  
  double operator()(const feature_type& feat, const double& grad) const
  {
    double& grad2 = const_cast<weight_set_type&>(grads2_).operator[](feat);
    const double eta = (grad2 == 0.0 ? eta0_ : eta0_ / utils::mathop::sqrt(grad2));
    
    grad2 += grad * grad;
    
    return eta;
  }
};

#endif
