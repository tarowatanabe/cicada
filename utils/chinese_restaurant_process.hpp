// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// Chinese Restaurant Process

#ifndef __UTILS__CHINESE_RESTAURANT_PROCESS__HPP__
#define __UTILS__CHINESE_RESTAURANT_PROCESS__HPP__ 1

#include <numeric>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <boost/functional/hash.hpp>

#include <utils/unordered_map.hpp>
#include <utils/slice_sampler.hpp>
#include <utils/mathop.hpp>

namespace utils
{

  template <typename Tp, typename Hash=boost::hash<Tp>, typename Pred=std::equal_to<Tp>, typename Alloc=std::allocator<Tp> >
  class chinese_restaurant_process
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Tp dish_type;
    
  public:
    chinese_restaurant_process()
      : table_size(),
	customer_size(),
	m_discount(0.9),
	m_strength(1.0),
	discount_prior_alpha(std::numeric_limits<double>::quiet_NaN()),
	discount_prior_beta(std::numeric_limits<double>::quiet_NaN()),
	strength_prior_shape(std::numeric_limits<double>::quiet_NaN()),
	strength_prior_rate(std::numeric_limits<double>::quiet_NaN())
    {}
    
    chinese_restaurant_process(const double& __discount,
			       const double& __strength)
      : table_size(),
	customer_size(),
	m_discount(__discount),
	m_strength(__strength),
	discount_prior_alpha(std::numeric_limits<double>::quiet_NaN()),
	discount_prior_beta(std::numeric_limits<double>::quiet_NaN()),
	strength_prior_shape(std::numeric_limits<double>::quiet_NaN()),
	strength_prior_rate(std::numeric_limits<double>::quiet_NaN())
    {}

    chinese_restaurant_process(const double& __discount_prior_alpha,
			       const double& __discount_prior_beta,
			       const double& __strength_prior_shape,
			       const double& __strength_prior_rate)
      : table_size(),
	customer_size(),
	m_discount(0.9),
	m_strength(1.0),
	discount_prior_alpha(__discount_prior_alpha),
	discount_prior_beta(__discount_prior_beta),
	strength_prior_shape(__strength_prior_shape),
	strength_prior_rate(__strength_prior_rate)
    {}
    
    chinese_restaurant_process(const double& __discount,
			       const double& __strength,
			       const double& __discount_prior_alpha,
			       const double& __discount_prior_beta,
			       const double& __strength_prior_shape,
			       const double& __strength_prior_rate)
      : table_size(),
	customer_size(),
	m_discount(__discount),
	m_strength(__strength),
	discount_prior_alpha(__discount_prior_alpha),
	discount_prior_beta(__discount_prior_beta),
	strength_prior_shape(__strength_prior_shape),
	strength_prior_rate(__strength_prior_rate)
    {}

  private:
    struct Location
    {
      typedef typename Alloc::template rebind<size_type>::other alloc_type;
      typedef std::vector<size_type, alloc_type> table_set_type;

      typedef typename table_set_type::const_iterator const_iterator;
      
      Location() : count(0) {}
      
      const_iterator begin() const { return tables.begin(); }
      const_iterator end() const { return tables.end(); }

      size_type size_customer() const { return count; }
      size_type size_table() const { return tables.size(); }
      
      size_type      count;
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

    bool has_discount_prior() const { return ! std::isnan(discount_prior_alpha); }
    bool has_strength_prior() const { return ! std::isnan(strength_prior_shape); }

    double& discount() { return m_discount; }
    double& strength() { return m_strength; }

    const double& discount() const { return m_discount; }
    const double& strength() const { return m_strength; }
    
    void clear()
    {
      table_size = 0;
      customer_size = 0;
      dishes.clear();
    }

    bool empty() const { return dishes.empty(); }
    
    size_type size_customer() const { return customer_size; }

    size_type size_table() const { return table_size; }
    
    size_type size_table(const dish_type& dish) const
    {
      typename dish_set_type::const_iterator diter = dishes.find(dish);
      
      return (diter == dishes.end() ? size_type(0) : diter->second.tables.size());
    }
    
    template <typename Sampler>
    bool increment(const dish_type& dish, const double& p0, Sampler& sampler)
    {
      location_type& loc = dishes[dish];
      
      bool shared = false;
      if (loc.count) {
	const double p_empty = (m_strength + table_size * m_discount) * p0;
	const double p_share = (loc.count - loc.tables.size() * m_discount);
	
	shared = sampler.select(p_empty, p_share);
      }
      
      if (shared) {
	double r = sampler.uniform() * (loc.count - loc.tables.size() * m_discount);
	
	typename location_type::table_set_type::iterator titer_end = loc.tables.end();
	for (typename location_type::table_set_type::iterator titer = loc.tables.begin(); titer != titer_end; ++ titer) {
	  r -= (*titer - m_discount);
	  
	  if (r <= 0.0) {
	    ++ (*titer);
	    break;
	  }
	}
      } else {
	loc.tables.push_back(1);
	++ table_size;
      }
      
      ++ loc.count;
      ++ customer_size;
      
      return ! shared;
    }
    
    template <typename Sampler>
    bool decrement(const dish_type& dish, Sampler& sampler)
    {
      typename dish_set_type::iterator diter = dishes.find(dish);
      
      if (diter == dishes.end())
	throw std::runtime_error("dish was not inserted?");
      
      location_type& loc = diter->second;
      
      if (loc.count == 1) {
	dishes.erase(diter);
	-- table_size;
	-- customer_size;
	return true;
      }
      
      bool erased = false;
      double r = sampler.uniform() * loc.count;
      -- loc.count;
      
      typename location_type::table_set_type::iterator titer_end = loc.tables.end();
      for (typename location_type::table_set_type::iterator titer = loc.tables.begin(); titer != titer_end; ++ titer) {
	r -= *titer;
	
	if (r <= 0.0) {
	  -- (*titer);
	  
	  if (! (*titer)) {
	    erased = true;
	    -- table_size;
	    loc.tables.erase(titer);
	  }
	  break;
	}
      }
      
      -- customer_size;
      return erased;
    }

    double prob(const dish_type& dish, const double& p0) const
    {
      typename dish_set_type::const_iterator diter = dishes.find(dish);
      
      const double r = table_size * m_discount + m_strength;
      
      if (diter == dishes.end())
	return r * p0 / (double(customer_size) + m_strength);
      else
	return (double(diter->second.count) - m_discount * diter->second.tables.size() + r * p0) / (double(customer_size) + m_strength);
    }

    double log_likelihood() const
    {
      return log_likelihood(m_discount, m_strength);
    }
    
    double log_likelihood(const double& discount, const double& strength) const
    {      
      double logprob = 0.0;
      
      if (has_discount_prior())
	logprob += utils::mathop::log_beta_density(discount, discount_prior_alpha, discount_prior_beta);
      
      if (has_strength_prior())
	logprob += utils::mathop::log_gamma_density(strength + discount, strength_prior_shape, strength_prior_rate);
      
      if (! customer_size) return logprob;
      
      if (discount > 0.0) {
	const double r = utils::mathop::lgamma(1.0 - discount);
	
	if (strength != 0.0)
          logprob += utils::mathop::lgamma(strength) - utils::mathop::lgamma(strength / discount);
	
	logprob += (- utils::mathop::lgamma(strength + customer_size)
		    + table_size * std::log(discount)
		    + utils::mathop::lgamma(strength / discount + table_size));
	
	typename dish_set_type::const_iterator diter_end = dishes.end();
	for (typename dish_set_type::const_iterator diter = dishes.begin(); diter != diter_end; ++ diter) {
	  const location_type& loc = diter->second;
	  
	  typename location_type::table_set_type::const_iterator titer_end = loc.tables.end();
	  for (typename location_type::table_set_type::const_iterator titer = loc.tables.begin(); titer != titer_end; ++ titer)
	    logprob += utils::mathop::lgamma(*titer - discount) - r;
	}
      } else if (discount == 0.0) {
	logprob += utils::mathop::lgamma(strength) + table_size * std::log(strength) - utils::mathop::lgamma(strength + table_size);
	
	typename dish_set_type::const_iterator diter_end = dishes.end();
	for (typename dish_set_type::const_iterator diter = dishes.begin(); diter != diter_end; ++ diter)
	  logprob += utils::mathop::lgamma(diter->second.tables.size());
      } else
	throw std::runtime_error("negative discount?");
      
      return logprob;
    }
    
    template <typename Sampler>
    void resample_hyperparameters(Sampler& sampler,
				  const size_type num_loop = 5,
				  const size_type num_iterations = 10)
    {
      if (! has_discount_prior() && ! has_strength_prior()) return;
      
      DiscountResampler discount_resampler(*this);
      StrengthResampler strength_resampler(*this);
      
      for (size_type iter = 0; iter < num_loop; ++iter) {
	if (has_strength_prior())
	  m_strength = slice_sampler(strength_resampler,
				     m_strength,
				     sampler,
				     - m_discount + std::numeric_limits<double>::min(),
				     std::numeric_limits<double>::infinity(),
				     0.0,
				     num_iterations,
				     100 * num_iterations);
	
	if (has_discount_prior()) 
	  m_discount = slice_sampler(discount_resampler,
				     m_discount,
				     sampler,
				     (m_strength < 0.0 ? - m_strength : 0.0) + std::numeric_limits<double>::min(),
				     1.0,
				     0.0,
				     num_iterations,
				     100 * num_iterations);
      }
      
      m_strength = slice_sampler(strength_resampler,
				 m_strength,
				 sampler,
				 - m_discount + std::numeric_limits<double>::min(),
				 std::numeric_limits<double>::infinity(),
				 0.0,
				 num_iterations,
				 100 * num_iterations);
    }

  private:    
    struct DiscountResampler
    {
      DiscountResampler(const chinese_restaurant_process& __crp) : crp(__crp) {}
      
      const chinese_restaurant_process& crp;
      
      double operator()(const double& proposed_discount) const
      {
	return crp.log_likelihood(proposed_discount, crp.strength());
      }
    };
    
    struct StrengthResampler
    {
      StrengthResampler(const chinese_restaurant_process& __crp) : crp(__crp) {}
      
      const chinese_restaurant_process& crp;
      
      double operator()(const double& proposed_strength) const
      {
	return crp.log_likelihood(crp.discount(), proposed_strength);
      }
    };
    
  private:
    size_type table_size;
    size_type customer_size;
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
