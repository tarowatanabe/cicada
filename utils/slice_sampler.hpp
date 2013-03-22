// -*- mode: c++ -*-
//
//  Copyright(C) 2012-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// taken from slice-sampler.h from py-cfg
//

//! slice-sampler.h is an MCMC slice sampler
//!
//! Mark Johnson, 30th September, 2012
//!
//! Based directly on the doubling procedure slice sampler described
//! in Neal (2003) "Slice Sampling", The Annals of Statistics 3:705-767.
//!
//! The primary user-callable functions are slice_sampler1d() and 
//! slice_sampler1dp(), which are defined at the bottom of this file.

#ifndef __UTILS__SLICE_SAMPLER__HPP__
#define __UTILS__SLICE_SAMPLER__HPP__ 1

#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#include <boost/lexical_cast.hpp>

namespace utils
{
  //! slice_sampler1d_type{} implements a 1-d slice sampler based on
  //! the procedure described in
  //! Neal (2003) "Slice Sampling", The Annals of Statistics 3:705-767.
  //
  template <typename F, typename Uniform01, typename Fn>
  struct slice_sampler1d_type {
    typedef size_t U;

    Uniform01& u01; //!< random number generator in range [0,1)
    const Fn& f;    //!< function to sample

    slice_sampler1d_type(Uniform01& u01, const Fn& f)
      : u01(u01), f(f)
    { }

    //! stepping_out() implements the "stepping out" procedure from Neal's Fig. 3.
    //
    void stepping_out(F x0, F y, F w, U m, F& l, F& r) {
      F u = u01();
      l = x0 - w*u;
      r = l + w;
      F v = u01();
      U j = m*v;
      U k = (m-1) * j;
      while (j > 0 && y < f(l)) {
	l -= w;
	j--;
      }
      while (k > 0 && y < f(r)) {
	r += w;
	k--;
      }
    }  // slice_sampler1d_type::stepping_out()

    //! doubling() implements the "doubling" procedure from Neal's Fig. 4.
    //
    void doubling(F x0, F y, F w, U p, F& l, F& r) {
      F u = u01();
      l = x0 - w*u;
      r = l + w;
      U k = p;
      while (k > 0 && (y < f(l) || y < f(r))) {
	F v = u01();
	if (v < 0.5)
	  l -= (r - l);
	else
	  r += (r - l);
	k--;
      }
    }  // slice_sampler1d_type::doubling()

    //! shrinkage() implements the "shrinkage" procedure from Neal's Fig. 5.
    //
    F shrinkage(F x0, F y, F w, F l, F r, bool always_accept) {
      // TRACE6(x0, y, w, l, r, always_accept);
      F lbar = l;
      F rbar = r;
      while (true) {
	F u = u01();
	F x1 = l + u*(rbar-lbar);
	F fx1 = f(x1);
	if (y < fx1 && (always_accept || acceptable(x0, x1, y, w, l, r)))
	  return x1;
	if (x1 < x0)
	  lbar = x1;
	else
	  rbar = x1;
	// TRACE5(x1, fx1, y, lbar, rbar);
      }
    } // slice_sampler1d_type::shrinkage()

    //! acceptable() implements the acceptance procedure from Neal's Fig. 6.
    //
    bool acceptable(F x0, F x1, F y, F w, F l, F r) const {
      bool d = false;
      while (r - l > 1.1*w) {
	//F m = (l+r)/2;
	F m = (l+r) * 0.5;
	if ((x0 < m && x1 >= m) || (x0 >= m && x1 < m))
	  d = true;
	if (x1 < m) 
	  r = m;
	else
	  l = m;
	if (d && y >= f(l) && y >= f(r)) {
	  //TRACE1(false);
	  return false;
	}
      }
      return true;
    }  // slice_sampler1d_type::acceptable()

    //! stepping_out_sample() computes one sample using the stepping-out
    //! procedure.
    //
    F stepping_out_sample(F x0, F w, U m) {
      F y = f(x0) + std::log(u01()+1e-100);
      F l, r;
      stepping_out(x0, y, w, m, l, r);
      F x1 = shrinkage(x0, y, w, l, r, true);
      return x1;
    } 

    //! doubling_sample() computes one sample using the doubling procedure.
    //
    F doubling_sample(F x0, F w, U p) {
      F y = f(x0) + std::log(u01()+1e-100);
      F l, r;
      doubling(x0, y, w, p, l, r);
      F x1 = shrinkage(x0, y, w, l, r, false);
      return x1;
    }

  }; // slice_sampler1d_type{}
  

  //! bounded_domain_function_type{}(x) returns f(x) if x is within the domain
  //!  bounds, and -infinity otherwise
  //
  template <typename F, typename Fn>
  struct bounded_domain_function_type {
    const Fn& f;      //!< original function
    F min_x, max_x;   //!< domain bounds
  
    bounded_domain_function_type(const Fn& f, F min_x, F max_x) 
      : f(f), min_x(min_x), max_x(max_x) { }
  
    F operator() (F x) const {
      if (min_x < x && x < max_x) {
	F fx = f(x);
	if (! std::isfinite(fx))
	  throw std::runtime_error(std::string("invalid function value:")
				   + " logF("+ boost::lexical_cast<std::string>(x) + ") = "
				   + boost::lexical_cast<std::string>(fx));
	return fx;
      }
      else
	return -std::numeric_limits<F>::infinity();      
    }
  };  // bounded_domain_function_type{}

  //! slice_sampler1d() implements the univariate slice sampler using the stepping-out procedure
  //! described in Neal (2003) "Slice Sampling", The Annals of Statistics 31(3), 705-767.
  //

  // renamed from slice_sampler1d to slice_sampler
  
  template <typename F, typename LogF, typename Uniform01>
  inline
  F slice_sampler(const LogF& logF,                                   //!< log of function to sample
		  F x0,                                               //!< starting point
		  Uniform01& u01,                                     //!< uniform [0,1) random number generator
		  F min_x = -std::numeric_limits<F>::infinity(),      //!< minimum value of support
		  F max_x = std::numeric_limits<F>::infinity(),       //!< maximum value of support
		  F w = 0.0,                                          //!< guess at initial width
		  size_t nsamples=1,                                   //!< number of samples to draw
		  size_t nsteps=32)                                    //!< number of width steps to permit
  {
    typedef size_t U;
    
    if (! std::isfinite(x0))
      throw std::runtime_error("invalid starting point: " + boost::lexical_cast<std::string>(x0));
    
    if (w <= 0.0) {                           // set w to a default width 
      if (min_x > -std::numeric_limits<F>::infinity() 
	  && max_x < std::numeric_limits<F>::infinity()) {
	//w = (max_x - min_x)/4;
	w = (max_x - min_x) * 0.25;
      }else {
	//w = std::max(((x0 < 0.0) ? -x0 : x0) /2, 1e-7);
	w = std::max(std::fabs(x0) * 0.5, 1e-7);
      }
    }

    if (! std::isfinite(w))
      throw std::runtime_error("invalid initial width: " + boost::lexical_cast<std::string>(w));
    
    typedef bounded_domain_function_type<F, LogF> BDLogF;
    BDLogF bdLogF(logF, min_x, max_x);
    slice_sampler1d_type<F, Uniform01, BDLogF> sampler(u01, bdLogF);
    
    for (U sample = 0; sample < nsamples; ++sample) {
      F x1 = sampler.stepping_out_sample(x0, w, nsteps);
      if (! std::isfinite(x1))
	throw std::runtime_error("invalid stepping out: " + boost::lexical_cast<std::string>(x1));
      
      w = std::fabs(1.5*(x1-x0));
      x0 = x1;
    }
    if (! std::isfinite(x0))
      throw std::runtime_error("invalid value: " + boost::lexical_cast<std::string>(x0));

    return x0;
  } // slice_sampler1d()

};

#endif
