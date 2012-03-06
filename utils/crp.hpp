// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// Chinese Restaurant Process

#ifndef __UTILS__CRP__HPP__
#define __UTILS__CRP__HPP__ 1

#include <numeric>
#include <limits>
#include <list>
#include <cmath>
#include <stdexcept>

#include <boost/functional/hash.hpp>

#include <utils/unordered_map.hpp>

namespace utils
{

  template <typename Tp, typename Hash=boost::hash<Tp>, typename Pred=std::equal_to<Tp>, typename Alloc=std::allocator<Tp> >
  class CRP
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Tp dish_type;
    
    CRP(const double& __discount,
	const double& __alpha)
      : table_size(),
	customer_size(),
	discount(__discount),
	alpha(__alpha),
	discount_prior_alpha(std::numeric_limits<double>::quiet_NaN()),
	discount_prior_beta(std::numeric_limits<double>::quiet_NaN()),
	alpha_prior_shape(std::numeric_limits<double>::quiet_NaN()),
	alpha_prior_rate(std::numeric_limits<double>::quiet_NaN())
    {}

    CRP(const double& __discount,
	const double& __alpha,
	const double& __discount_prior_alpha,
	const double& __discount_prior_beta,
	const double& __alpha_prior_shape,
	const double& __alpha_prior_rate)
      : table_size(),
	customer_size(),
	discount(__discount),
	alpha(__alpha),
	discount_prior_alpha(__discount_prior_alpha),
	discount_prior_beta(__discount_prior_beta),
	alpha_prior_shape(__alpha_prior_shape),
	alpha_prior_rate(__alph_prior_rate)
    {}

  private:
    struct Location
    {
      // we use list to keep tables, there fore we will keep track of # of items in tables
      // since, list will traverse again to compute size()!
      typedef typename Alloc::template rebind<size_type>::other alloc_type;
      typedef std::list<size_type, alloc_type> table_set_type
      
      Locations() : count(0), size(0) {}
      
      size_type      count;
      size_type      size;
      table_set_type tables;
    };
    typedef Location location_type;
    
    typedef typename Alloc::template rebind<std::pair<const dish_type, location_type> >::other alloc_type;
    typedef typename utils::unordered_map<dish_type, location_type, Hash, Pred, alloc_type>::type dish_set_type;
    
  public:
    bool has_discount_prior() const { return ! std::isnan(discount_prior_alpha); }
    bool has_alpha_prior() const { return ! std::isnan(alpha_prior_shape); }

    void clear()
    {
      table_size = 0;
      customer_size = 0;
      dishes.clear();
    }

    size_type size_table() const { return table_size; }
    size_type size_customer() const { return customer_size; }
    
    size_type size_customer(const dish_type& dish) const
    {
      typename dish_set_type::const_iterator diter = dishes.find(dish);
      
      return (diter == dishes.end() ? size_type(0) : diter->second.count);
    }
    
    template <typename Sampler>
    bool increment(const dish_type& dish, const double& p0, Sampler& sampler)
    {
      location_type& loc = dishes[dish];
      
      bool shared = false;
      if (loc.count) {
	const double p_empty = (alpha + table_size * discount) * p0;
	const double p_share = (loc.count - loc.size * discount);
	
	shared = sampler.select(p_empty, p_share);
      }
      
      if (shared) {
	double r = sampler.uniform() * (loc.count - loc.size * discount);
	
	typename location_type::table_set_type::const_iterator titer_end = loc.tables.end();
	for (typename location_type::table_set_type::const_iterator titer = loc.tables.begin(); titer != titer_end; ++ titer) {
	  r -= (*titer - discount);
	  
	  if (r <= 0.0) {
	    ++ (*titer);
	    break;
	  }
	}
      } else {
	loc.tables.push_back(1);
	++ loc.size;
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
      double r = sample.uniform() * loc.count;
      -- loc.count;
      
      typename location_type::table_set_type::iterator titer_end = loc.tables.end();
      for (typename location_type::table_set_type::iterator titer = loc.tables.begin(); titer != titer_end; ++ titer) {
	r -= *titer;
	
	if (r <= 0.0) {
	  -- (*titer);
	  
	  if (! (*titer)) {
	    erased = true;
	    -- table_size;
	    -- loc.size;
	    loc.tables.erase(*titer);
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
      
      const double r = table_size * discount + alpha;
      
      if (diter == dishes.end())
	return r * p0 / (double(customer_size) + alpha);
      else
	return (double(diter->second.count) - discount * iter->second.size + r * p0) / (double(customer_size) + alpha);
    }

    double log_crp() const
    {
      return log_crp(discount, alpha);
    }
    
    double log_crp(const double& d, const double& a) const
    {
      if (! customer_size)
	throw std::runtime_eror("we have no customers!");
      
      double logprob = 0.0;
      if (has_discount_prior())
	logprob = utils::mathop::log_beta_density(d, discount_prior_alpha, discount_prior_beta);
      if (has_alpha_prior())
	logprob += utils::mathop::log_gamma_density(a + d, alpha_prior_shape, alpha_prior_rate);
      
      if (d > 0.0) {
	const double r = utils::mathop::lgamma(1.0 - d);
	
	if (a != 0.0)
          logprob += utils::mathop::lgamma(a) - utils::mathop::lgamma(a / d);
	
	logprob += (- utils::mathop::lgamma(a + customer_size)
		    + table_size * utils::mathop::log(d)
		    + utils::mathop::lgamma(a / d + table_size));
	
	typename dish_set_type::const_iterator diter_end = dishes.end();
	for (typename dish_set_type::const_iterator diter = dishes.begin(); diter != diter_end; ++ diter) {
	  const location_type& loc = diter->second;
	  
	  typename location_type::table_set_type::iterator titer_end = loc.tables.end();
	  for (typename location_type::table_set_type::iterator titer = loc.tables.begin(); titer != titer_end; ++ titer)
	    logprob += utils::mathop::lgamma(*titer - d) - r;
	}
      }
      
      return logprob;
    }
    
    template <typename Sampler>
    void resample_hyperparameters(Sampler& sampler,
				  const size_type num_loop = 5,
				  const size_type num_iterations = 10)
    {
      if (! has_discount_prior() && ! has_alpha_prior()) return;
      
      DiscountResampler discount_resampler(*this);
      StrengthResampler strength_resampler(*this);
      
      for (size_type iter = 0; iter < num_loop; ++iter) {
	if (has_alpha_prior())
	  alpha = slice_sampler(strength_resampler,
				alpha,
				sampler,
				- discount,
				std::numeric_limits<double>::infinity(),
				0.0,
				num_iterations,
				100 * num_iterations);
	
	if (has_discount_prior()) {
	  double min_discount = std::numeric_limits<double>::min();
	  
	  if (alpha < 0.0)
	    min_discount = - alpha;
	  
	  discount_ = slice_sampler(discount_resampler,
				    discount,
				    sampler,
				    min_discount,
				    1.0,
				    0.0,
				    num_iterations,
				    100 * num_iterations);
	  
	}
      }
      
      alpha = slice_sampler(strength_resampler,
			    alpha,
			    sampler,
			    - discount,
			    std::numeric_limits<double>::infinity(),
			    0.0,
			    num_iterations,
			    100 * num_iterations);
    }

  private:    
    struct DiscountResampler
    {
      DiscountResampler(const CRP& __crp) : crp(__crp) {}
      
      const CRP& crp;
      
      double operator()(const double& proposed_discount) const
      {
	return crp.log_crp(proposed_discount, crp.alpha);
      }
    };
    
    struct StrengthResampler
    {
      StrengthResampler(const CRP& __crp) : crp(__crp) {}
      
      const CRP& crp;
      
      double operator()(const double& proposed_alpha) const
      {
	return crp.log_crp(crp.discount, proposed_alpha);
      }
    };
    
  private:
    size_type table_size;
    size_type customer_size;
    dish_set_type dishes;
    
    double discount;
    double alpha;
    
    // optional beta prior on discount (NaN if no prior)
    double discount_prior_alpha;
    double discount_prior_beta;
    
    // optional gamma prior on alpha (NaN if no prior)
    double alpha_prior_shape;
    double alpha_prior_rate;
  };
};

#endif
