// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//


#ifndef __UTILS__RESTAURANT__HPP__
#define __UTILS__RESTAURANT__HPP__ 1

#include <numeric>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>

#include <utils/unordered_map.hpp>
#include <utils/slice_sampler.hpp>
#include <utils/mathop.hpp>
#include <utils/pyp_parameter.hpp>

// Chinese Restaurant Process
//
// inspired by restauraht.hh
//

//
// Chinese Restaurant Process with optional customer and table tracking
// Copyright Trevor Cohn 2009, University of Edinburgh
// Table tracking algorithm courtesy of Mark Johnson, and is described in:
//
//      A Note on the Implementation of Hierarchical Dirichlet Processes
//      Phil Blunsom, Trevor Cohn, Sharon Goldwater and Mark Johnson
//      ACL 2009 (short paper)
//

namespace utils
{

  template <typename Tp, typename Hash=boost::hash<Tp>, typename Pred=std::equal_to<Tp>, typename Alloc=std::allocator<Tp> >
  class restaurant
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Tp dish_type;
    typedef pyp_parameter parameter_type;
    
  public:
    restaurant()
      : tables(),
	customers(),
	parameter() {}
    
    restaurant(const parameter_type& __parameter)
      : tables(),
	customers(),
	parameter(__parameter)
    {
      parameter.verify_parameters();
    }
    
    restaurant(const double& __discount,
	       const double& __strength)
      : tables(),
	customers(),
	parameter(__discount, __strength)
    {
      parameter.verify_parameters();
    }

    restaurant(const double& __discount_prior_alpha,
	       const double& __discount_prior_beta,
	       const double& __strength_prior_shape,
	       const double& __strength_prior_rate)
      : tables(),
	customers(),
	parameter(__discount_prior_alpha,
		  __discount_prior_beta,
		  __strength_prior_shape,
		  __strength_prior_rate)
    {
      parameter.verify_parameters();
    }
    
    restaurant(const double& __discount,
	       const double& __strength,
	       const double& __discount_prior_alpha,
	       const double& __discount_prior_beta,
	       const double& __strength_prior_shape,
	       const double& __strength_prior_rate)
      : tables(),
	customers(),
	parameter(__discount,
		  __strength,
		  __discount_prior_alpha,
		  __discount_prior_beta,
		  __strength_prior_shape,
		  __strength_prior_rate)
    {
      parameter.verify_parameters();
    }

  private:
    struct Location
    {
      typedef typename Alloc::template rebind<size_type>::other alloc_type;
      typedef std::vector<size_type, alloc_type> table_set_type;

      typedef typename table_set_type::const_iterator const_iterator;
      
      Location() : customers(0) {}
      
      const_iterator begin() const { return tables.begin(); }
      const_iterator end() const { return tables.end(); }

      size_type size_customer() const { return customers; }
      size_type size_table() const { return tables.size(); }
      
      size_type      customers;
      table_set_type tables;
    };
    typedef Location location_type;
    
    typedef typename Alloc::template rebind<std::pair<const dish_type, location_type> >::other alloc_type;
    typedef typename utils::unordered_map<dish_type, location_type, Hash, Pred, alloc_type>::type dish_set_type;

  public:
    typedef typename dish_set_type::key_type       key_type;
    typedef typename dish_set_type::mapped_type    mapped_type;
    typedef typename dish_set_type::value_type     value_type;
    typedef typename dish_set_type::const_iterator const_iterator;
    typedef typename dish_set_type::const_iterator iterator;
    
  public:
    const_iterator begin() const { return dishes.begin(); }
    const_iterator end() const { return dishes.end(); }

    bool has_discount_prior() const { return parameter.has_discount_prior(); }
    bool has_strength_prior() const { return parameter.has_strength_prior(); }

    double& discount() { return parameter.discount; }
    double& strength() { return parameter.strength; }

    const double& discount() const { return parameter.discount; }
    const double& strength() const { return parameter.strength; }
    
    void clear()
    {
      tables = 0;
      customers = 0;
      dishes.clear();
    }

    bool empty() const { return dishes.empty(); }
    
    size_type size_customer() const { return customers; }

    size_type size_table() const { return tables; }
    
    size_type size_table(const dish_type& dish) const
    {
      typename dish_set_type::const_iterator diter = dishes.find(dish);
      
      return (diter == dishes.end() ? size_type(0) : diter->second.tables.size());
    }

    size_type size_customer(const dish_type& dish) const
    {
      typename dish_set_type::const_iterator diter = dishes.find(dish);
      
      return (diter == dishes.end() ? size_type(0) : diter->second.customers);
    }

    void swap(restaurant& x)
    {
      std::swap(tables, x.tables);
      std::swap(customers, x.customers);
      dishes.swap(x.dishes);
      parameter.swap(x.parameter);
    }

    
    template <typename Sampler>
    bool increment(const dish_type& dish, const double& p0, Sampler& sampler, const double temperature=1.0)
    {
      location_type& loc = dishes[dish];
      
      bool existing = false;
      if (loc.customers) {
	if (temperature == 1.0) {
	  const double p_base = (parameter.strength + tables * parameter.discount) * p0;
	  const double p_gen  = (loc.customers - loc.tables.size() * parameter.discount);
	  
	  existing = sampler.bernoulli(p_gen / (p_base + p_gen));
	} else {
	  const double p_base = std::pow((parameter.strength + tables * parameter.discount) * p0, 1.0 / temperature);
	  const double p_gen  = std::pow((loc.customers - loc.tables.size() * parameter.discount), 1.0 / temperature);
	  
	  existing = sampler.bernoulli(p_gen / (p_base + p_gen));
	}
      }
      
      if (existing) {
	double r = sampler.uniform() * (loc.customers - loc.tables.size() * parameter.discount);
	
	typename location_type::table_set_type::iterator titer_end = loc.tables.end();
	for (typename location_type::table_set_type::iterator titer = loc.tables.begin(); titer != titer_end; ++ titer) {
	  r -= (*titer - parameter.discount);
	  
	  if (r <= 0.0) {
	    ++ (*titer);
	    break;
	  }
	}
      } else {
	loc.tables.push_back(1);
	++ tables;
      }
      
      ++ loc.customers;
      ++ customers;
      
      return ! existing;
    }
    
    template <typename Sampler>
    bool decrement(const dish_type& dish, Sampler& sampler)
    {
      typename dish_set_type::iterator diter = dishes.find(dish);
      
      if (diter == dishes.end())
	throw std::runtime_error("dish was not inserted?");
      
      location_type& loc = diter->second;
      
      if (loc.customers == 1) {
	dishes.erase(diter);
	-- tables;
	-- customers;
	return true;
      }
      
      bool erased = false;
      double r = sampler.uniform() * loc.customers;
      -- loc.customers;
      
      typename location_type::table_set_type::iterator titer_end = loc.tables.end();
      for (typename location_type::table_set_type::iterator titer = loc.tables.begin(); titer != titer_end; ++ titer) {
	r -= *titer;
	
	if (r <= 0.0) {
	  -- (*titer);
	  
	  if (! (*titer)) {
	    erased = true;
	    -- tables;
	    loc.tables.erase(titer);
	  }
	  break;
	}
      }
      
      -- customers;
      return erased;
    }
    
    template <typename P>
    P prob(const dish_type& dish, const P& p0) const
    {
      typename dish_set_type::const_iterator diter = dishes.find(dish);
      
      if (diter == dishes.end())
	return P(tables * parameter.discount + parameter.strength) * p0 / P(customers + parameter.strength);
      else
	return (P(diter->second.customers - parameter.discount * diter->second.tables.size()) + P(tables * parameter.discount + parameter.strength) * p0) / P(customers + parameter.strength);
    }
    
    double log_likelihood() const
    {
      return log_likelihood(parameter.discount, parameter.strength);
    }
    
    // http://en.wikipedia.org/wiki/Chinese_restaurant_process
    double log_likelihood(const double& discount, const double& strength) const
    {      
      double logprob = 0.0;
      
      if (has_discount_prior())
	logprob += utils::mathop::log_beta_density(discount, parameter.discount_alpha, parameter.discount_beta);
      
      if (has_strength_prior())
	logprob += utils::mathop::log_gamma_density(strength + discount, parameter.strength_shape, parameter.strength_rate);
      
      if (! customers) return logprob;
      
      if (discount > 0.0) {
	if (strength == 0.0)
	  logprob += tables * std::log(discount) + utils::mathop::lgamma(tables) - utils::mathop::lgamma(customers);
	else {
	  logprob += utils::mathop::lgamma(strength) - utils::mathop::lgamma(strength + customers);
	  logprob += tables * std::log(discount) + utils::mathop::lgamma(strength / discount + tables) - utils::mathop::lgamma(strength / discount);
	}
	
	const double lg = utils::mathop::lgamma(1.0 - discount);
	
	typename dish_set_type::const_iterator diter_end = dishes.end();
	for (typename dish_set_type::const_iterator diter = dishes.begin(); diter != diter_end; ++ diter) {
	  const location_type& loc = diter->second;
	  
	  typename location_type::table_set_type::const_iterator titer_end = loc.tables.end();
	  for (typename location_type::table_set_type::const_iterator titer = loc.tables.begin(); titer != titer_end; ++ titer)
	    logprob += utils::mathop::lgamma(*titer - discount) - lg;
	}
      } else if (discount == 0.0) {
	logprob += utils::mathop::lgamma(strength) + tables * std::log(strength) - utils::mathop::lgamma(strength + tables);
	
	typename dish_set_type::const_iterator diter_end = dishes.end();
	for (typename dish_set_type::const_iterator diter = dishes.begin(); diter != diter_end; ++ diter)
	  logprob += utils::mathop::lgamma(diter->second.tables.size());
      } else
	throw std::runtime_error("negative discount?");
      
      return logprob;
    }

  private:    
    struct DiscountSampler
    {
      DiscountSampler(const restaurant& __crp) : crp(__crp) {}
      
      const restaurant& crp;
      
      double operator()(const double& proposed_discount) const
      {
	return crp.log_likelihood(proposed_discount, crp.strength());
      }
    };
    
    struct StrengthSampler
    {
      StrengthSampler(const restaurant& __crp) : crp(__crp) {}
      
      const restaurant& crp;
      
      double operator()(const double& proposed_strength) const
      {
	return crp.log_likelihood(crp.discount(), proposed_strength);
      }
    };

  public:
    template <typename Sampler>
    double sample_log_x(Sampler& sampler, const double& discount, const double& strength) const
    {
      return std::log(sample_x(sampler, discount, strength));
    }
    
    template <typename Sampler>
    double sample_x(Sampler& sampler, const double& discount, const double& strength) const
    {
      return (customers > 1 ? sampler.beta(strength + 1, customers - 1) : 1.0);
    }

    template <typename Sampler>
    double sample_y(Sampler& sampler, const double& discount, const double& strength) const
    {
      size_type y = 0;
      
      for (size_type i = 1; i < tables; ++ i)
	y += sampler.bernoulli(strength / (strength + discount * i));
      
      return y;
    }

    template <typename Sampler>
    double sample_y_inv(Sampler& sampler, const double& discount, const double& strength) const
    {
      size_type y = 0;
      
      for (size_type i = 1; i < tables; ++ i)
	y += 1 - sampler.bernoulli(strength / (strength + discount * i));
      
      return y;
    }
    
    template <typename Sampler>
    double sample_z_inv(Sampler& sampler, const double& discount, const double& strength) const
    {
      size_type z = 0;
      
      typename dish_set_type::const_iterator diter_end = dishes.end();
      for (typename dish_set_type::const_iterator diter = dishes.begin(); diter != diter_end; ++ diter) {
	const location_type& loc = diter->second;
	
	typename location_type::table_set_type::const_iterator titer_end = loc.tables.end();
	for (typename location_type::table_set_type::const_iterator titer = loc.tables.begin(); titer != titer_end; ++ titer)
	  for (size_type j = 1; j < *titer; ++ j)
	    z += 1 - sampler.bernoulli(double(j - 1) / (j - discount));
      }
      
      return z;
    }

    template <typename Sampler>
    double sample_strength(Sampler& sampler, const double& discount, const double& strength) const
    {
      if (! has_strength_prior())
	throw std::runtime_error("no strength prior");
      
      const double x = sample_log_x(sampler, discount, strength);
      const double y = sample_y(sampler, discount, strength);
      
      return sampler.gamma(parameter.strength_shape + y, parameter.strength_rate - x);
    }
    
    template <typename Sampler>
    double sample_discount(Sampler& sampler, const double& discount, const double& strength) const
    {
      if (! has_discount_prior())
	throw std::runtime_error("no discount prior");
      
      const double y = sample_y_inv(sampler, discount, strength);
      const double z = sample_z_inv(sampler, discount, strength);
      
      return sampler.beta(parameter.discount_alpha + y, parameter.discount_beta + z);
    }
    
    template <typename Sampler>
    double expectation_strength(Sampler& sampler, const double& discount, const double& strength) const
    {
      const double x = sample_log_x(sampler, discount, strength);
      const double y = sample_y(sampler, discount, strength);
      
      if (has_strength_prior())
	return (parameter.strength_shape + y) / (parameter.strength_rate - x);
      else
	return - y / x;
    }
    
    template <typename Sampler>
    double expectation_discount(Sampler& sampler, const double& discount, const double& strength) const
    {
      const double y = sample_y_inv(sampler, discount, strength);
      const double z = sample_z_inv(sampler, discount, strength);
      
      if (has_discount_prior()) {
	const double a = parameter.discount_alpha + y;
	const double b = parameter.discount_beta  + z;
	
	return a / (a + b);
      } else
	return y / (y + z);
    }
    
    bool verify_parameters()
    {
      return parameter.verify_parameters();
    }
    
    template <typename Sampler>
    void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
    {
      if (! has_discount_prior() && ! has_strength_prior()) return;
      
      for (int iter = 0; iter != num_loop; ++ iter) {
	if (has_strength_prior())
	  parameter.strength = sample_strength(sampler, parameter.discount, parameter.strength);
	
	if (has_discount_prior()) 
	  parameter.discount = sample_discount(sampler, parameter.discount, parameter.strength);
      }
      
      if (has_strength_prior())
	parameter.strength = sample_strength(sampler, parameter.discount, parameter.strength);
    }
    
    template <typename Sampler>
    void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
    {
      if (! has_discount_prior() && ! has_strength_prior()) return;
      
      DiscountSampler discount_sampler(*this);
      StrengthSampler strength_sampler(*this);

      for (int iter = 0; iter != num_loop; ++ iter) {
	if (has_strength_prior())
	  parameter.strength = slice_sampler(strength_sampler,
					     parameter.strength,
					     sampler,
					     - parameter.discount + std::numeric_limits<double>::min(),
					     std::numeric_limits<double>::infinity(),
					     0.0,
					     num_iterations,
					     100 * num_iterations);
	
	if (has_discount_prior()) 
	  parameter.discount = slice_sampler(discount_sampler,
					     parameter.discount,
					     sampler,
					     (parameter.strength < 0.0 ? - parameter.strength : 0.0) + std::numeric_limits<double>::min(),
					     1.0,
					     0.0,
					     num_iterations,
					     100 * num_iterations);
      }
      
      if (has_strength_prior())
	parameter.strength = slice_sampler(strength_sampler,
					   parameter.strength,
					   sampler,
					   - parameter.discount + std::numeric_limits<double>::min(),
					   std::numeric_limits<double>::infinity(),
					   0.0,
					   num_iterations,
					   100 * num_iterations);
    }
    
  private:
    size_type      tables;
    size_type      customers;
    dish_set_type  dishes;
    parameter_type parameter;
  };
};

namespace std
{
  template <typename T, typename H, typename P, typename A>
  inline
  void swap(utils::restaurant<T,H,P,A>& x, utils::restaurant<T,H,P,A>& y)
  {
    x.swap(y);
  }
};

#endif
