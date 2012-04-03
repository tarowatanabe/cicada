// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//


#ifndef __UTILS__RESTAURANT_FLOOR__HPP__
#define __UTILS__RESTAURANT_FLOOR__HPP__ 1

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

// Chinese Restaurant_Floor Process
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

// multiple-floor CRP
//
// @inproceedings{WooTeh2009a,
// Author =        {F. Wood and Y. W. Teh},
// Title =         {A Hierarchical Nonparametric {B}ayesian Approach 
//                 to Statistical Language Model Domain Adaptation},
// Booktitle =     {Proceedings of the International Conference on 
//                 Artificial Intelligence and Statistics},
// Volume =        {12},
// Year =          {2009}}

namespace utils
{

  template <typename Tp, typename Hash=boost::hash<Tp>, typename Pred=std::equal_to<Tp>, typename Alloc=std::allocator<Tp> >
  class restaurant_floor
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Tp dish_type;
    
  public:
    restaurant_floor()
      : tables(),
	customers(),
	m_discount(0.9),
	m_strength(1.0),
	discount_prior_alpha(std::numeric_limits<double>::quiet_NaN()),
	discount_prior_beta(std::numeric_limits<double>::quiet_NaN()),
	strength_prior_shape(std::numeric_limits<double>::quiet_NaN()),
	strength_prior_rate(std::numeric_limits<double>::quiet_NaN())
    {
      verify_parameters();
    }
    
    restaurant_floor(const double& __discount,
		     const double& __strength)
      : tables(),
	customers(),
	m_discount(__discount),
	m_strength(__strength),
	discount_prior_alpha(std::numeric_limits<double>::quiet_NaN()),
	discount_prior_beta(std::numeric_limits<double>::quiet_NaN()),
	strength_prior_shape(std::numeric_limits<double>::quiet_NaN()),
	strength_prior_rate(std::numeric_limits<double>::quiet_NaN())
    {
      verify_parameters();
    }

    restaurant_floor(const double& __discount_prior_alpha,
		     const double& __discount_prior_beta,
		     const double& __strength_prior_shape,
		     const double& __strength_prior_rate)
      : tables(),
	customers(),
	m_discount(0.9),
	m_strength(1.0),
	discount_prior_alpha(__discount_prior_alpha),
	discount_prior_beta(__discount_prior_beta),
	strength_prior_shape(__strength_prior_shape),
	strength_prior_rate(__strength_prior_rate)
    {
      verify_parameters();
    }
    
    restaurant_floor(const double& __discount,
		     const double& __strength,
		     const double& __discount_prior_alpha,
		     const double& __discount_prior_beta,
		     const double& __strength_prior_shape,
		     const double& __strength_prior_rate)
      : tables(),
	customers(),
	m_discount(__discount),
	m_strength(__strength),
	discount_prior_alpha(__discount_prior_alpha),
	discount_prior_beta(__discount_prior_beta),
	strength_prior_shape(__strength_prior_shape),
	strength_prior_rate(__strength_prior_rate)
    {
      verify_parameters();
    }

  private:
    struct Location
    {
      typedef size_type count_type;
      typedef size_type floor_type;
      typedef std::pair<size_type, size_type> value_type;
      typedef typename Alloc::template rebind<value_type>::other alloc_type;
      
      typedef std::vector<value_type, alloc_type> table_set_type;

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

    bool has_discount_prior() const { return ! std::isnan(discount_prior_alpha) && ! std::isnan(discount_prior_beta); }
    bool has_strength_prior() const { return ! std::isnan(strength_prior_shape) && ! std::isnan(strength_prior_rate); }

    double& discount() { return m_discount; }
    double& strength() { return m_strength; }

    const double& discount() const { return m_discount; }
    const double& strength() const { return m_strength; }
    
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
    
    template <typename PriorIterator, typename LambdaIterator, typename Sampler>
    std::pair<size_type, bool> increment(const dish_type& dish,
					 PriorIterator first,
					 PriorIterator last,
					 LambdaIterator lambda,
					 Sampler& sampler,
					 const double temperature=1.0)
    {
      location_type& loc = dishes[dish];
      
      const double p0 = std::inner_product(first, last, lambda, 0.0);
      size_type floor = 0;
      
      bool existing = false;
      if (loc.customers) {
	if (temperature == 1.0) {
	  const double p_base = (m_strength + tables * m_discount) * p0;
	  const double p_gen  = (loc.customers - loc.tables.size() * m_discount);
	  
	  existing = sampler.bernoulli(p_gen / (p_base + p_gen));
	} else {
	  const double p_base = std::pow((m_strength + tables * m_discount) * p0, 1.0 / temperature);
	  const double p_gen  = std::pow((loc.customers - loc.tables.size() * m_discount), 1.0 / temperature);
	  
	  existing = sampler.bernoulli(p_gen / (p_base + p_gen));
	}
      }
      
      if (existing) {
	double r = sampler.uniform() * (loc.customers - loc.tables.size() * m_discount);
	
	typename location_type::table_set_type::iterator titer_end = loc.tables.end();
	for (typename location_type::table_set_type::iterator titer = loc.tables.begin(); titer != titer_end; ++ titer) {
	  r -= titer->first - m_discount;
	  
	  if (r <= 0.0) {
	    ++ titer->first;
	    floor = titer->second;
	    break;
	  }
	}
      } else {
	// sample what floor...
	const size_type num_floor = std::distance(first, last);

	if (num_floor > 1) {
	  double r = sampler.uniform() * p0;
	  for (/**/; floor != num_floor; ++ floor, ++ first, ++ lambda) {
	    r -= (*first) * (*lambda);
	    
	    if (r <= 0.0) break;
	  }
	}
	
	loc.tables.push_back(std::make_pair(1, floor));
	++ tables;
      }
      
      ++ loc.customers;
      ++ customers;
      
      return std::make_pair(floor, ! existing);
    }
    
    template <typename Sampler>
    std::pair<size_type, bool> decrement(const dish_type& dish, Sampler& sampler)
    {
      typename dish_set_type::iterator diter = dishes.find(dish);
      
      if (diter == dishes.end())
	throw std::runtime_error("dish was not inserted?");
      
      location_type& loc = diter->second;
      
      if (loc.customers == 1) {
	const size_type floor = loc.tables.front().second;
	dishes.erase(diter);
	-- tables;
	-- customers;
	return std::make_pair(floor, true);
      }
      
      size_type floor = 0;
      bool erased = false;
      double r = sampler.uniform() * loc.customers;
      -- loc.customers;
      
      typename location_type::table_set_type::iterator titer_end = loc.tables.end();
      for (typename location_type::table_set_type::iterator titer = loc.tables.begin(); titer != titer_end; ++ titer) {
	r -= titer->first;
	
	if (r <= 0.0) {
	  floor = titer->second;
	  -- titer->first;
	  
	  if (! titer->first) {
	    erased = true;
	    -- tables;
	    loc.tables.erase(titer);
	  }
	  break;
	}
      }
      
      -- customers;
      return std::make_pair(floor, erased);
    }
    
    template <typename PriorIterator, typename LambdaIterator>
    typename std::iterator_traits<PriorIterator>::value_type prob(const dish_type& dish, PriorIterator first, PriorIterator last, LambdaIterator lambda) const
    {
      typedef typename std::iterator_traits<PriorIterator>::value_type P;
      typename dish_set_type::const_iterator diter = dishes.find(dish);
      
      const P p0 = std::inner_product(first, last, lambda, P(0));
      
      if (diter == dishes.end())
	return P(tables * m_discount + m_strength) * p0 / P(customers + m_strength);
      else
	return (P(diter->second.customers - m_discount * diter->second.tables.size()) + P(tables * m_discount + m_strength) * p0) / P(customers + m_strength);
    }
    
    double log_likelihood() const
    {
      return log_likelihood(m_discount, m_strength);
    }
    
    // http://en.wikipedia.org/wiki/Chinese_restaurant_process
    double log_likelihood(const double& discount, const double& strength) const
    {      
      double logprob = 0.0;
      
      if (has_discount_prior())
	logprob += utils::mathop::log_beta_density(discount, discount_prior_alpha, discount_prior_beta);
      
      if (has_strength_prior())
	logprob += utils::mathop::log_gamma_density(strength + discount, strength_prior_shape, strength_prior_rate);
      
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
	    logprob += utils::mathop::lgamma(titer->first - discount) - lg;
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
      DiscountSampler(const restaurant_floor& __crp) : crp(__crp) {}
      
      const restaurant_floor& crp;
      
      double operator()(const double& proposed_discount) const
      {
	return crp.log_likelihood(proposed_discount, crp.strength());
      }
    };
    
    struct StrengthSampler
    {
      StrengthSampler(const restaurant_floor& __crp) : crp(__crp) {}
      
      const restaurant_floor& crp;
      
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
      
      return sampler.gamma(strength_prior_shape + y, strength_prior_rate - x);
    }
    
    template <typename Sampler>
    double sample_discount(Sampler& sampler, const double& discount, const double& strength) const
    {
      if (! has_discount_prior())
	throw std::runtime_error("no discount prior");
      
      const double y = sample_y_inv(sampler, discount, strength);
      const double z = sample_z_inv(sampler, discount, strength);
      
      return sampler.beta(discount_prior_alpha + y, discount_prior_beta + z);
    }
    
    template <typename Sampler>
    double expectation_strength(Sampler& sampler, const double& discount, const double& strength) const
    {
      const double x = sample_log_x(sampler, discount, strength);
      const double y = sample_y(sampler, discount, strength);
      
      if (has_strength_prior())
	return (strength_prior_shape + y) / (strength_prior_rate - x);
      else
	return - y / x;
    }
    
    template <typename Sampler>
    double expectation_discount(Sampler& sampler, const double& discount, const double& strength) const
    {
      const double y = sample_y_inv(sampler, discount, strength);
      const double z = sample_z_inv(sampler, discount, strength);
      
      if (has_discount_prior()) {
	const double a = discount_prior_alpha + y;
	const double b = discount_prior_beta  + z;
	
	return a / (a + b);
      } else
	return y / (y + z);
    }
    
    bool verify_parameters()
    {
      if (m_discount < 0.0 || m_discount >= 1.0)
	throw std::runtime_error("invalid discount: " + boost::lexical_cast<std::string>(m_discount));
      
      if (m_strength <= - m_discount)
	throw std::runtime_error("invalid strength: " + boost::lexical_cast<std::string>(m_strength));
      
      return true;
    }
    
    template <typename Sampler>
    void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
    {
      if (! has_discount_prior() && ! has_strength_prior()) return;
      
      for (int iter = 0; iter != num_loop; ++ iter) {
	if (has_strength_prior())
	  m_strength = sample_strength(sampler, m_discount, m_strength);
	
	if (has_discount_prior()) 
	  m_discount = sample_discount(sampler, m_discount, m_strength);
      }
      
      if (has_strength_prior())
	m_strength = sample_strength(sampler, m_discount, m_strength);
    }
    
    template <typename Sampler>
    void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
    {
      if (! has_discount_prior() && ! has_strength_prior()) return;
      
      DiscountSampler discount_sampler(*this);
      StrengthSampler strength_sampler(*this);

      for (int iter = 0; iter != num_loop; ++ iter) {
	if (has_strength_prior())
	  m_strength = slice_sampler(strength_sampler,
				     m_strength,
				     sampler,
				     - m_discount + std::numeric_limits<double>::min(),
				     std::numeric_limits<double>::infinity(),
				     0.0,
				     num_iterations,
				     100 * num_iterations);
	
	if (has_discount_prior()) 
	  m_discount = slice_sampler(discount_sampler,
				     m_discount,
				     sampler,
				     (m_strength < 0.0 ? - m_strength : 0.0) + std::numeric_limits<double>::min(),
				     1.0,
				     0.0,
				     num_iterations,
				     100 * num_iterations);
      }
      
      if (has_strength_prior())
	m_strength = slice_sampler(strength_sampler,
				   m_strength,
				   sampler,
				   - m_discount + std::numeric_limits<double>::min(),
				   std::numeric_limits<double>::infinity(),
				   0.0,
				   num_iterations,
				   100 * num_iterations);
    }
    
  private:
    size_type tables;
    size_type customers;
    dish_set_type dishes;
    
    double m_discount;
    double m_strength;
    
    // optional beta prior on discount (NaN if no prior)
    double discount_prior_alpha;
    double discount_prior_beta;
    
    // optional gamma prior on strength (NaN if no prior)
    double strength_prior_shape;
    double strength_prior_rate;
  };
};

#endif
