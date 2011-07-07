// -*- mode: c++ -*-

#ifndef __LINEARPP__TRON__HPP__
#define __LINEARPP__TRON__HPP__ 1

#include <algorithm>

#ifdef __cplusplus
extern "C" {
#endif

  extern double dnrm2_(int *, double *, int *);
  extern double ddot_(int *, double *, int *, double *, int *);
  extern int daxpy_(int *, double *, double *, int *, double *, int *);
  extern int dscal_(int *, double *, double *, int *);

#ifdef __cplusplus
}
#endif

namespace linearpp
{
  class Tron
  {
  public:
    Tron(const double __eps=0.1, const int __max_iter = 1000)
      : eps(__eps), max_iter(__max_iter) {}

    template <typename Function, typename Gradient, typename Hv hv>
    double tron(const size_t num_features, Function func, Gradient grad, Hv hv, double* w)
    {      
      // Parameters for updating the iterates.
      double eta0 = 1e-4, eta1 = 0.25, eta2 = 0.75;

      // Parameters for updating the trust region size delta.
      double sigma1 = 0.25, sigma2 = 0.5, sigma3 = 4;

      int n = num_features;
      int i, cg_iter;
      double delta, snorm, one=1.0;
      double alpha, f, fnew, prered, actred, gs;
      int search = 1, iter = 1, inc = 1;
      double *s = new double[n];
      double *r = new double[n];
      double *w_new = new double[n];
      double *g = new double[n];
      
      for (i=0; i<n; i++)
	w[i] = 0;

      f = fun(w);
      grad(w, g);
      delta = dnrm2_(&n, g, &inc);
      double gnorm1 = delta;
      double gnorm = gnorm1;

      if (gnorm <= eps*gnorm1)
	search = 0;

      iter = 1;

      while (iter <= max_iter && search)
	{
	  cg_iter = trcg(num_features, hv, delta, g, s, r);
	  
	  memcpy(w_new, w, sizeof(double)*n);
	  daxpy_(&n, &one, s, &inc, w_new, &inc);

	  gs = ddot_(&n, g, &inc, s, &inc);
	  prered = -0.5*(gs-ddot_(&n, s, &inc, r, &inc));
	  fnew = fun(w_new);

	  // Compute the actual reduction.
	  actred = f - fnew;

	  // On the first iteration, adjust the initial step bound.
	  snorm = dnrm2_(&n, s, &inc);
	  if (iter == 1)
	    delta = std::min(delta, snorm);

	  // Compute prediction alpha*snorm of the step.
	  if (fnew - f - gs <= 0)
	    alpha = sigma3;
	  else
	    alpha = std::max(sigma1, -0.5*(gs/(fnew - f - gs)));

	  // Update the trust region bound according to the ratio of actual to predicted reduction.
	  if (actred < eta0*prered)
	    delta = std::min(std::max(alpha, sigma1)*snorm, sigma2*delta);
	  else if (actred < eta1*prered)
	    delta = std::max(sigma1*delta, std::min(alpha*snorm, sigma2*delta));
	  else if (actred < eta2*prered)
	    delta = std::max(sigma1*delta, std::min(alpha*snorm, sigma3*delta));
	  else
	    delta = std::max(delta, std::min(alpha*snorm, sigma3*delta));

	  info("iter %2d act %5.3e pre %5.3e delta %5.3e f %5.3e |g| %5.3e CG %3d\n", iter, actred, prered, delta, f, gnorm, cg_iter);

	  if (actred > eta0*prered)
	    {
	      iter++;
	      memcpy(w, w_new, sizeof(double)*n);
	      f = fnew;
	      grad(w, g);

	      gnorm = dnrm2_(&n, g, &inc);
	      if (gnorm <= eps*gnorm1)
		break;
	    }
	  if (f < -1.0e+32)
	    {
	      info("warning: f < -1.0e+32\n");
	      break;
	    }
	  if (fabs(actred) <= 0 && prered <= 0)
	    {
	      info("warning: actred and prered <= 0\n");
	      break;
	    }
	  if (fabs(actred) <= 1.0e-12*fabs(f) &&
	      fabs(prered) <= 1.0e-12*fabs(f))
	    {
	      info("warning: actred and prered too small\n");
	      break;
	    }
	}

      delete[] g;
      delete[] r;
      delete[] w_new;
      delete[] s;

      return f;
    }

    template <typename Hv>
    int trcg(const size_t num_features, Hv hv, double delta, double* g, double* s, double* r)
    {
      int i, inc = 1;
      int n = num_features;
      double one = 1;
      double *d = new double[n];
      double *Hd = new double[n];
      double rTr, rnewTrnew, alpha, beta, cgtol;

      for (i=0; i<n; i++)
	{
	  s[i] = 0;
	  r[i] = -g[i];
	  d[i] = r[i];
	}
      cgtol = 0.1*dnrm2_(&n, g, &inc);

      int cg_iter = 0;
      rTr = ddot_(&n, r, &inc, r, &inc);
      while (1)
	{
	  if (dnrm2_(&n, r, &inc) <= cgtol)
	    break;
	  cg_iter++;
	  hv(d, Hd);

	  alpha = rTr/ddot_(&n, d, &inc, Hd, &inc);
	  daxpy_(&n, &alpha, d, &inc, s, &inc);
	  if (dnrm2_(&n, s, &inc) > delta)
	    {
	      info("cg reaches trust region boundary\n");
	      alpha = -alpha;
	      daxpy_(&n, &alpha, d, &inc, s, &inc);

	      double std = ddot_(&n, s, &inc, d, &inc);
	      double sts = ddot_(&n, s, &inc, s, &inc);
	      double dtd = ddot_(&n, d, &inc, d, &inc);
	      double dsq = delta*delta;
	      double rad = sqrt(std*std + dtd*(dsq-sts));
	      if (std >= 0)
		alpha = (dsq - sts)/(std + rad);
	      else
		alpha = (rad - std)/dtd;
	      daxpy_(&n, &alpha, d, &inc, s, &inc);
	      alpha = -alpha;
	      daxpy_(&n, &alpha, Hd, &inc, r, &inc);
	      break;
	    }
	  alpha = -alpha;
	  daxpy_(&n, &alpha, Hd, &inc, r, &inc);
	  rnewTrnew = ddot_(&n, r, &inc, r, &inc);
	  beta = rnewTrnew/rTr;
	  dscal_(&n, &beta, d, &inc);
	  daxpy_(&n, &one, r, &inc, d, &inc);
	  rTr = rnewTrnew;
	}

      delete[] d;
      delete[] Hd;

      return(cg_iter);
    }
    
  private:
    double eps;
    int    max_iter;
  };
};

#endif
