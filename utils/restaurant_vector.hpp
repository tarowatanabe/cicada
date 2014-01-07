// -*- mode: c++ -*-
//
//  Copyright(C) 2012-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//


#ifndef __UTILS__RESTAURANT_VECTOR__HPP__
#define __UTILS__RESTAURANT_VECTOR__HPP__ 1

#include <numeric>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <boost/lexical_cast.hpp>

#include <utils/slice_sampler.hpp>
#include <utils/mathop.hpp>
#include <utils/pyp_parameter.hpp>
#include <utils/table_count.hpp>

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
  template <typename Alloc=std::allocator<size_t> >
  class restaurant_vector
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef size_type dish_type;
    typedef pyp_parameter parameter_type;
    
  public:
    restaurant_vector()
      : tables(),
	customers(),
	parameter() {}
    
    restaurant_vector(const parameter_type& __parameter)
      : tables(),
	customers(),
	parameter(__parameter)
    {
      parameter.verify_parameters();
    }
    
    restaurant_vector(const double& __discount,
		      const double& __strength)
      : tables(),
	customers(),
	parameter(__discount, __strength)
    {
      parameter.verify_parameters();
    }

    restaurant_vector(const double& __discount_alpha,
		      const double& __discount_beta,
		      const double& __strength_shape,
		      const double& __strength_rate)
      : tables(),
	customers(),
	parameter(__discount_alpha,
		  __discount_beta,
		  __strength_shape,
		  __strength_rate)
    {
      parameter.verify_parameters();
    }
    
    restaurant_vector(const double& __discount,
		      const double& __strength,
		      const double& __discount_alpha,
		      const double& __discount_beta,
		      const double& __strength_shape,
		      const double& __strength_rate)
      : tables(),
	customers(),
	parameter(__discount,
		  __strength,
		  __discount_alpha,
		  __discount_beta,
		  __strength_shape,
		  __strength_rate)
    {
      parameter.verify_parameters();
    }

    
  private:
    struct Location
    {
      typedef typename Alloc::template rebind<size_type>::other alloc_type;

      typedef utils::table_count<size_type, 1, alloc_type> count_set_type;
      
      typedef typename count_set_type::histogram_type::const_iterator const_iterator;
      
      Location() : counts() {}
      
      const_iterator begin() const { return counts[0].begin(); }
      const_iterator end() const { return counts[0].end(); }

      size_type size_customer() const { return counts.customers(); }
      size_type size_table() const { return counts.tables(); }
      
      bool empty() const { return counts.customers() == 0; }
      
      void clear()
      {
	counts.clear();
      }

      void swap(Location& x)
      {
	counts.swap(x.counts);
      }
      
      count_set_type counts;
    };
    typedef Location location_type;
    
    typedef typename Alloc::template rebind<location_type >::other alloc_type;
    typedef std::vector<location_type, alloc_type> dish_set_type;

  public:
    typedef typename dish_set_type::value_type     value_type;
    typedef typename dish_set_type::const_iterator const_iterator;
    typedef typename dish_set_type::const_iterator iterator;
    
  public:
    const_iterator begin() const { return dishes.begin(); }
    const_iterator end() const { return dishes.end(); }
    
    location_type& operator[](size_type x) { return dishes[x]; }
    const location_type& operator[](size_type x) const { return dishes[x]; }

    void reserve(size_type x) { dishes.reserve(x); }
    void resize(size_type x) { dishes.resize(x); }

    bool has_discount_prior() const { return parameter.has_discount_prior(); }
    bool has_strength_prior() const { return parameter.has_strength_prior(); }

    double& discount() { return parameter.discount; }
    double& strength() { return parameter.strength; }

    const double& discount() const { return parameter.discount; }
    const double& strength() const { return parameter.strength; }

    const parameter_type& parameters() const { return parameter; }
    
    void clear()
    {
      tables = 0;
      customers = 0;
      dishes.clear();
    }

    bool empty() const { return dishes.empty(); }
    size_type size() const { return dishes.size(); }
    
    size_type size_customer() const { return customers; }
    size_type size_table() const { return tables; }
    
    size_type size_table(const dish_type dish) const
    {
      return (dish < dishes.size() ? dishes[dish].size_table() : size_type(0));
    }

    size_type size_customer(const dish_type dish) const
    {
      return (dish < dishes.size() ? dishes[dish].size_customer() : size_type(0));
    }

    void swap(restaurant_vector& x)
    {
      std::swap(tables, x.tables);
      std::swap(customers, x.customers);
      dishes.swap(x.dishes);
      parameter.swap(x.parameter);
    }

    template <typename Mapping>
    void permute(const Mapping& mapping)
    {
      dish_set_type dishes_new(mapping.size());
      
      for (size_type i = 0; i != mapping.size(); ++ i) 
	if (mapping[i] < dishes.size())
	  dishes_new[i].swap(dishes[mapping[i]]);

      for (size_type i = 0; i != dishes.size(); ++ i)
	if (! dishes[i].empty())
	  throw std::runtime_error("non-empty dish i=" + boost::lexical_cast<std::string>(i));
      
      while (! dishes_new.empty() && dishes_new.back().empty())
	dishes_new.pop_back();
      
      dishes.swap(dishes_new);
    }

    template <typename Sampler>
    bool increment_existing(const dish_type& dish, Sampler& sampler)
    {
      if (dish >= dishes.size())
	dishes.resize(dish + 1);
      
      location_type& loc = dishes[dish];
      
      if (loc.empty()) {
	++ customers;
	++ tables;
	loc.counts.increment_new();
	return true;
      } else {
	++ customers;
	loc.counts.increment_existing(parameter.discount, sampler);
	return false;
      }
    }
    
    template <typename Sampler>
    bool increment_new(const dish_type& dish, Sampler& sampler)
    {
      if (dish >= dishes.size())
	dishes.resize(dish + 1);
      
      ++ customers;
      ++ tables;
      dishes[dish].counts.increment_new();
      return true;
    }
    
    template <typename Sampler>
    bool increment(const dish_type dish, const double& p0, Sampler& sampler, const double temperature=1.0)
    {
      if (dish >= dishes.size())
	dishes.resize(dish + 1);
      
      location_type& loc = dishes[dish];
      
      bool existing = false;
      if (loc.size_customer()) {
	if (temperature == 1.0) {
	  const double p_base = (parameter.strength + tables * parameter.discount) * p0;
	  const double p_gen  = (loc.size_customer() - loc.size_table() * parameter.discount);
	  
	  existing = sampler.bernoulli(p_gen / (p_base + p_gen));
	} else {
	  const double p_base = std::pow((parameter.strength + tables * parameter.discount) * p0, 1.0 / temperature);
	  const double p_gen  = std::pow((loc.size_customer() - loc.size_table() * parameter.discount), 1.0 / temperature);
	  
	  existing = sampler.bernoulli(p_gen / (p_base + p_gen));
	}
      }
      
      if (! existing) {
	++ customers;
	++ tables;
	loc.counts.increment_new();
	return true;
      } else {
	++ customers;
	loc.counts.increment_existing(parameter.discount, sampler);
	return false;
      }
    }
    
    template <typename Sampler>
    bool decrement(const dish_type dish, Sampler& sampler)
    {
      if (dish >= dishes.size())
	throw std::runtime_error("restaurant_vector: dish was not inserted?");
      
      location_type& loc = dishes[dish];

      const bool removed = loc.counts.decrement(sampler).second;
      
      -- customers;
      tables -= removed;
      
      return removed;
    }
    
    template <typename P>
    P prob(const dish_type dish, const P& p0) const
    {
      if (dish >= dishes.size())
	return P(tables * parameter.discount + parameter.strength) * p0 / P(customers + parameter.strength);
      else
	return (P(dishes[dish].size_customer() - parameter.discount * dishes[dish].size_table()) + P(tables * parameter.discount + parameter.strength) * p0) / P(customers + parameter.strength);
    }
    
    template <typename P>
    P prob(const P& p0) const
    {
      return P(tables * parameter.discount + parameter.strength) * p0 / P(customers + parameter.strength);
    }
    
    template <typename P>
    std::pair<P, bool> prob_model(const dish_type dish, const P& p0) const
    {
      if (dish >= dishes.size() || dishes[dish].empty())
	return std::make_pair(P(tables * parameter.discount + parameter.strength) * p0 / P(customers + parameter.strength), false);
      else
	return std::make_pair((P(dishes[dish].size_customer() - parameter.discount * dishes[dish].size_table()) + P(tables * parameter.discount + parameter.strength) * p0) / P(customers + parameter.strength), true);
    }
    
    double log_likelihood() const
    {
      return log_likelihood(parameter.discount, parameter.strength);
    }
    
    // http://en.wikipedia.org/wiki/Chinese_restaurant_process
    double log_likelihood(const double& discount, const double& strength) const
    {      
      double logprob = parameter.log_likelihood();
      
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
	for (typename dish_set_type::const_iterator diter = dishes.begin(); diter != diter_end; ++ diter) 
	  if (! diter->empty()) {
	    const location_type& loc = *diter;
	    
	    typename location_type::const_iterator titer_end = loc.end();
	    for (typename location_type::const_iterator titer = loc.begin(); titer != titer_end; ++ titer)
	      logprob += (utils::mathop::lgamma(titer->first - discount) - lg) * titer->second;
	  }
      } else if (discount == 0.0) {
	logprob += utils::mathop::lgamma(strength) + tables * std::log(strength) - utils::mathop::lgamma(strength + tables);
	
	typename dish_set_type::const_iterator diter_end = dishes.end();
	for (typename dish_set_type::const_iterator diter = dishes.begin(); diter != diter_end; ++ diter)
	  if (! diter->empty())
	    logprob += utils::mathop::lgamma(diter->size_table());
      } else
	throw std::runtime_error("negative discount?");
      
      return logprob;
    }
    
    void prune()
    {
      dish_set_type(dishes).swap(dishes);
    }

  private:    
    struct DiscountSampler
    {
      DiscountSampler(const restaurant_vector& __crp) : crp(__crp) {}
      
      const restaurant_vector& crp;
      
      double operator()(const double& proposed_discount) const
      {
	return crp.log_likelihood(proposed_discount, crp.strength());
      }
    };
    
    struct StrengthSampler
    {
      StrengthSampler(const restaurant_vector& __crp) : crp(__crp) {}
      
      const restaurant_vector& crp;
      
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
	const location_type& loc = *diter;
	
	typename location_type::const_iterator titer_end = loc.end();
	for (typename location_type::const_iterator titer = loc.begin(); titer != titer_end; ++ titer)
	  for (size_type i = 0; i != titer->second; ++ i)
	    for (size_type j = 1; j < titer->first; ++ j)
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
					     32 * num_iterations);
	
	if (has_discount_prior()) 
	  parameter.discount = slice_sampler(discount_sampler,
					     parameter.discount,
					     sampler,
					     (parameter.strength < 0.0 ? - parameter.strength : 0.0) + std::numeric_limits<double>::min(),
					     1.0,
					     0.0,
					     num_iterations,
					     32 * num_iterations);
      }
      
      if (has_strength_prior())
	parameter.strength = slice_sampler(strength_sampler,
					   parameter.strength,
					   sampler,
					   - parameter.discount + std::numeric_limits<double>::min(),
					   std::numeric_limits<double>::infinity(),
					   0.0,
					   num_iterations,
					   32 * num_iterations);
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
  template <typename A>
  inline
  void swap(utils::restaurant_vector<A>& x, utils::restaurant_vector<A>& y)
  {
    x.swap(y);
  }
  
};


#endif
