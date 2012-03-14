// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MATHOP__HPP__
#define __UTILS__MATHOP__HPP__ 1

#ifdef HAVE_TR1_CMATH
#include <tr1/cmath>
#endif

#include <cmath>
#include <cfloat>

#include <algorithm>

#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include <utils/config.hpp>

namespace utils
{
  namespace mathop
  {
    
    template <typename Tp>
    inline
    Tp factorial(unsigned n)
    {
      using namespace boost::math::policies;
      typedef policy<domain_error<errno_on_error>,
	pole_error<errno_on_error>,
	overflow_error<errno_on_error>,
	rounding_error<errno_on_error>,
	evaluation_error<errno_on_error>
	> policy_type;
      
      const Tp ret = boost::math::factorial<Tp>(n, policy_type());
      return (std::isfinite(ret) ? ret : boost::numeric::bounds<Tp>::highest());
    }
    
    template <typename Tp>
    inline
    Tp log(Tp value)
    {
      const Tp ret = std::log(value);
      return (std::isfinite(ret) ? ret : boost::numeric::bounds<Tp>::lowest());
    }
    
    template <typename Tp>
    inline
    Tp exp(Tp value)
    {
      const Tp ret = std::exp(value);
      return (std::isfinite(ret) ? ret : boost::numeric::bounds<Tp>::highest());
    }
    
    template <typename Tp>
    inline
    Tp pow(Tp base, const Tp x)
    {
      const Tp ret = std::pow(base, x);
      return (std::isfinite(ret) ? ret : boost::numeric::bounds<Tp>::highest());
    }

    template <typename Tp>
    inline
    Tp sqrt(Tp value)
    {
      const Tp ret = std::sqrt(value);
      return (std::isfinite(ret) ? ret : 0.0);
    }
    
    template <typename Tp>
    inline
    Tp logsum(Tp x, Tp y)
    {
#ifdef HAVE_TR1_CMATH
      if (x <= boost::numeric::bounds<Tp>::lowest())
        return y;
      else if (y <= boost::numeric::bounds<Tp>::lowest())
        return x;
      else
	return std::max(x, y) + std::tr1::log1p(mathop::exp(std::min(x, y) - std::max(x, y)));
#else
      using namespace boost::math::policies;
      typedef policy<domain_error<errno_on_error>,
	pole_error<errno_on_error>,
	overflow_error<errno_on_error>,
	rounding_error<errno_on_error>,
	evaluation_error<errno_on_error>
	> policy_type;
    
      if (x <= boost::numeric::bounds<Tp>::lowest())
        return y;
      else if (y <= boost::numeric::bounds<Tp>::lowest())
        return x;
      else
	return std::max(x, y) + boost::math::log1p(mathop::exp(std::min(x, y) - std::max(x, y)), policy_type());
#endif
    }
    
    template <typename Tp>
    inline
    Tp digamma(Tp x)
    {
      // this works as if log(value) is performed...
      using namespace boost::math::policies;
      typedef policy<domain_error<errno_on_error>,
	pole_error<errno_on_error>,
	overflow_error<errno_on_error>,
	rounding_error<errno_on_error>,
	evaluation_error<errno_on_error>
	> policy_type;
    
      const Tp ret = boost::math::digamma(x, policy_type());
      return (std::isfinite(ret) ? ret : boost::numeric::bounds<Tp>::lowest());
    }

    template <typename Tp>
    inline
    Tp lgamma(Tp x)
    {
#ifdef HAVE_TR1_CMATH
      return std::tr1::lgamma(x);
#else
      return boost::math::lgamma(x);
#endif
    }
    
    template <typename Tp>
    inline
    Tp log_poisson(size_t x, const Tp& lambda) 
    {
      return std::log(lambda) * x - lgamma(x + 1) - lambda;
    }

    template <typename Tp>
    inline
    Tp log_gamma(const Tp& x) 
    {
      return lgamma(x);
    }
    
    template <typename Tp>
    inline
    Tp log_beta(const Tp& x, const Tp& y) 
    {
      return lgamma(x) + lgamma(y) - lgamma(x + y);
    }
    
    template <typename Tp>
    inline
    Tp log_gamma_density(const Tp& x, const Tp& shape, const Tp& rate) 
    {
      return (shape - Tp(1)) * std::log(x) - shape * std::log(rate) - x / rate - lgamma(shape);
    }
    
    template <typename Tp>
    inline
    Tp log_beta_density(const Tp& x, const Tp& alpha, const Tp& beta) 
    {
      return (alpha - Tp(1)) * std::log(x) + (beta - Tp(1)) * std::log(Tp(1) - x) - log_beta(alpha, beta);
    }

    
  };
};

#endif
