// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// taken from slice_sampler.h from py-cfg
// 
//! slice-sampler.h is an MCMC slice sampler
//!
//! Mark Johnson, 1st August 2008 

#ifndef __UTILS__SLICE_SAMPLER__HPP__
#define __UTILS__SLICE_SAMPLER__HPP__ 1

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

namespace utils
{
  //! slice_sampler_rfc_type{} returns the value of a user-specified
  //! function if the argument is within range, or - infinity otherwise
  //

  namespace impl
  {
    template <typename F, typename Fn, typename U>
    struct slice_sampler_rfc_type {
      F min_x;
      F max_x;
      const Fn& f;
      U max_nfeval;
      U nfeval;
      
      slice_sampler_rfc_type(F min_x, F max_x, const Fn& f, U max_nfeval) 
	: min_x(min_x), max_x(max_x), f(f), max_nfeval(max_nfeval), nfeval(0) { }
      
      F operator() (const F& x) {
	if (min_x < x && x < max_x) {
	  
	  ++ nfeval;
	  if (nfeval > max_nfeval)
	    throw std::runtime_error("exceed max function evaluation");
	  
	  const F fx = f(x);
	  
	  if (! std::isfinite(fx))
	    throw std::runtime_error("invlaid function value");
	  
	  return fx;
	} else
	  return - std::numeric_limits<F>::infinity();
      }
    };
  };  // slice_sampler_rfc_type{}
  
  //! slice_sampler1d() implements the univariate "range doubling" slice sampler
  //! described in Neal (2003) "Slice Sampling", The Annals of Statistics 31(3), 705-767.
  //
  
  template <typename F, typename LogF, typename Uniform01>
  F slice_sampler(const LogF& logF0,               //!< log of function to sample
		  F x,                             //!< starting point
		  Uniform01& u01,                  //!< uniform [0,1) random number generator
		  F min_x = -std::numeric_limits<F>::infinity(),  //!< minimum value of support
		  F max_x = std::numeric_limits<F>::infinity(),   //!< maximum value of support
		  F w = 0.0,                       //!< guess at initial width
		  size_t nsamples=1,             //!< number of samples to draw
		  size_t max_nfeval=200)         //!< max number of function evaluations
  {
    typedef size_t U;
    impl::slice_sampler_rfc_type<F,LogF,U> logF(min_x, max_x, logF0, max_nfeval);

    if (! std::isfinite(x))
      throw std::runtime_error("x should be finite");
    
    if (w <= 0.0) {                           // set w to a default width 
      if (min_x > - std::numeric_limits<F>::infinity() && max_x < std::numeric_limits<F>::infinity())
	w = (max_x - min_x) / 4;
      else
	w = std::max(std::fabs(x) / 4, 0.1);
    }
    
    if (! std::isfinite(w))
      throw std::runtime_error("invalid w");
    
    F logFx = logF(x);
    for (U sample = 0; sample < nsamples; ++ sample) {
      F logY = logFx + std::log(u01() + 1e-100);     //! slice logFx at this value
      
      if (! std::isfinite(logY))
	throw std::runtime_error("invalid logY=" + boost::lexical_cast<std::string>(logY) + " logFx=" + boost::lexical_cast<std::string>(logFx));
      
      F xl = x - w * u01();                   //! lower bound on slice interval
      F logFxl = logF(xl);
      F xr = xl + w;                          //! upper bound on slice interval
      F logFxr = logF(xr);

      while (logY < logFxl || logY < logFxr)  // doubling procedure
	if (u01() < 0.5) 
	  logFxl = logF(xl -= xr - xl);
	else
	  logFxr = logF(xr += xr - xl);
	
      F xl1 = xl;
      F xr1 = xr;
      for (;;) {                          // shrinking procedure
	F x1 = xl1 + u01() * (xr1 - xl1);
	if (logY < logF(x1)) {
	  F xl2 = xl;                         // acceptance procedure
	  F xr2 = xr; 
	  bool d = false;
	  while (xr2 - xl2 > 1.1 * w) {
	    F xm = (xl2 + xr2) / 2;
	    if ((x < xm && x1 >= xm) || (x >= xm && x1 < xm))
	      d = true;
	    if (x1 < xm)
	      xr2 = xm;
	    else
	      xl2 = xm;
	    if (d && logY >= logF(xl2) && logY >= logF(xr2))
	      goto unacceptable;
	  }
	  x = x1;
	  goto acceptable;
	}
	goto acceptable;
      unacceptable:
	if (x1 < x)                           // rest of shrinking procedure
	  xl1 = x1;
	else 
	  xr1 = x1;
      }
    acceptable:
      w = (4 * w + (xr1 - xl1)) / 5;              // update width estimate
    }
    return x;
  }
  
};

#endif
