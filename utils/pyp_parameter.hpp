// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__PYP_PARAMETER__HPP__
#define __UTILS__PYP_PARAMETER__HPP__ 1

// a structrue which hold pyp parameters...

#include <stdexcept>
#include <limits>

#include <utils/mathop.hpp>

namespace utils
{
  struct pyp_parameter
  {
    pyp_parameter()
      : discount(0.0),
	strength(1.0),
	discount_alpha(std::numeric_limits<double>::quiet_NaN()),
	discount_beta(std::numeric_limits<double>::quiet_NaN()),
	strength_shape(std::numeric_limits<double>::quiet_NaN()),
	strength_rate(std::numeric_limits<double>::quiet_NaN()) {}
    
    pyp_parameter(const double& __discount,
		  const double& __strength)
      : discount(__discount),
	strength(__strength),
	discount_alpha(std::numeric_limits<double>::quiet_NaN()),
	discount_beta(std::numeric_limits<double>::quiet_NaN()),
	strength_shape(std::numeric_limits<double>::quiet_NaN()),
	strength_rate(std::numeric_limits<double>::quiet_NaN())
    {
      verify_parameters();
    }

    pyp_parameter(const double& __discount_alpha,
		  const double& __discount_beta,
		  const double& __strength_shape,
		  const double& __strength_rate)
      : discount(0.9),
	strength(1.0),
	discount_alpha(__discount_alpha),
	discount_beta(__discount_beta),
	strength_shape(__strength_shape),
	strength_rate(__strength_rate)
    {
      verify_parameters();
    }
    
    pyp_parameter(const double& __discount,
		  const double& __strength,
		  const double& __discount_alpha,
		  const double& __discount_beta,
		  const double& __strength_shape,
		  const double& __strength_rate)
      : discount(__discount),
	strength(__strength),
	discount_alpha(__discount_alpha),
	discount_beta(__discount_beta),
	strength_shape(__strength_shape),
	strength_rate(__strength_rate)
    {
      verify_parameters();
    }   
    
    bool has_discount_prior() const { return ! std::isnan(discount_alpha) && ! std::isnan(discount_beta); }
    bool has_strength_prior() const { return ! std::isnan(strength_shape) && ! std::isnan(strength_rate); }
    
    bool verify_parameters()
    {
      if (discount < 0.0 || discount >= 1.0)
	throw std::runtime_error("invalid discount: " + boost::lexical_cast<std::string>(discount));
      
      if (strength <= - discount)
	throw std::runtime_error("invalid strength: " + boost::lexical_cast<std::string>(strength));
      
      return true;
    }
    
    double log_likelihood() const
    {
      return log_likelihood(discount, strength);
    }
    
    double log_likelihood(const double& discount, const double& strength) const
    {
      double logprob = 0.0;
      
      if (has_discount_prior())
	logprob += utils::mathop::log_beta_density(discount, discount_alpha, discount_beta);
      
      if (has_strength_prior())
	logprob += utils::mathop::log_gamma_density(strength + discount, strength_shape, strength_rate);
      
      return logprob;
    }

    double discount;
    double strength;
    
    double discount_alpha;
    double discount_beta;
    double strength_shape;
    double strength_rate;
  };
};

#endif
