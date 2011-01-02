// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MATHOP__HPP__
#define __UTILS__MATHOP__HPP__ 1

#include <cmath>
#include <cfloat>

#include <algorithm>

#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/numeric/conversion/bounds.hpp>

namespace utils
{
  namespace mathop
  {
    
    template <typename Tp>
    inline
    double factorial(Tp n)
    {
      double ret = 1.0;
      for (/**/; n > 0; -- n)
	ret *= double(n);
      return (std::isfinite(ret) ? ret : boost::numeric::bounds<double>::highest());
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
  };
};

#endif
