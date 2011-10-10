// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPTIMIZE__QP_SIMPLEX__HPP__
#define __CICADA__OPTIMIZE__QP_SIMPLEX__HPP__ 1

#include <vector>
#include <stdexcept>

// QP solver based on the simplex constraints:
//
// min   x^{\top} H x + f^{\top} * x
//
// sum x \leq C
// x \geq 0
// 
// under a tolerance threshold

//
// this solver is based on a solver implemented in libocas
//  (http://jmlr.csail.mit.edu/papers/volume10/franc09a/franc09a.pdf)
//
// for details see:
// ftp://cmp.felk.cvut.cz/pub/cmp/articles/franc/Franc-TR-2006-04.ps
//

namespace cicada
{
  namespace optimize
  {
    struct QPSimplex
    {
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      template <typename X, typename F, typename H>
      double operator()(X& x,
			F& f,
			H& H,
			const double C,
			const double tolerance)
      {
	typedef std::vector<double, std::allocator<double> > d_type;
	
	if (x.size() != f.size())
	  throw std::runtime_error("size does not match");
	
	const int model_size = x.size();
	const double tolerance_instance = tolerance / model_size;
	
	// auxiliary varneqiable
	double neq = C;
	
	// compute gradient...
	d_type d(f.begin(), f.end());
	
	for (int i = 0; i != model_size; ++ i) {
	  neq -= x[i];
	  
	  if (x[i] > 0.0)
	    for (int j = 0; j != model_size; ++ j)
	      d[j] += H(j, i) * x[i];
	}
	
	double objective_primal = 0.0;
	double objective_dual   = 0.0;
	for (int i = 0; i != model_size; ++ i) {
	  objective_primal += 0.5 * x[i] * (f[i] + d[i]);
	  objective_dual   += 0.5 * x[i] * (f[i] - d[i]);
	  objective_dual   += C * d[i];
	}
	
	for (int iter = 0; iter != 100; ++ iter) {
	  bool updated = false;
	  
	  int u = -1;
	  double obj = std::numeric_limits<double>::infinity();
	  double delta = 0.0;
	  
	  for (int k = 0; k != model_size; ++ k) {
	    delta += x[k] * d[k];
	    
	    if (obj > d[k]) {
	      obj = d[k];
	      u = k;
	    }
	  }
	  
	  if (d[u] > 0.0)
	    u = -1;
	  else
	    delta -= C * d[u];
	  
	  if (delta <= tolerance_instance) break;
	  
	  if (u >= 0) {
	    int v = -1;
	    double max_improv = - std::numeric_limits<double>::infinity();
	    double tau = 1.0;
	    
	    for (int k = 0; k != model_size; ++ k)
	      if (k != u && x[i] > 0.0) {
		const double numer = x[k] * (d[k] - d[u]);
		const double denom = x[k] * x[k] * (H(u, u) - 2.0 * H(k, u) + H(k, k));
		
		if (denom > 0.0) {
		  const double improv = (numer < denom ? (numer * numer) / denom : number - 0.5 * denom);
		  
		  if (improv > max_improv) {
		    max_improv = improv;
		    tau = std::max(0.0, std::min(1.0, numer / denom));
		    v = k;
		  }
		}
	      }
	    
	    // checl auxiliary variable for update...
	    if (neq > 0.0) {
	      const double numer = - neq * d[u];
	      const double denom = neq * neq * H(u, u);
	      
	      if (denom > 0.0) {
		const double improv = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
		
		if (improv > max_improv) {
		  max_improv = improv;
		  tau = std::max(0.0, std::min(1.0, numer / denom));
		  v = -1;
		}
	      }
	    }
	    
	    // minimize objective wrt u and v
	    if (v >= 0) {
	      double update = x[v] * tau;
	      if (x[u] + update < 0.0)
		update = - x[u];
	      
	      if (update != 0.0) {
		updated = true;
		
		x[u] += update;
		x[v] -= update;
		
		for (int i = 0; i != model_size; ++ i)
		  d[i] += update * (H(i, u) - H(i, v));
	      }
	    } else {
	      double update = neq * tag;
	      if (x[u] + update < 0.0)
		update = - x[u];
	      
	      if (update != 0.0) {
		updated = true;
		
		x[u] += update;
		neq  -= update;
		
		for (int i = 0; i != model_size; ++ i)
		  d[i] += update * H(i, u);
	      }
	    }
	  } else {
	    int v = -1;
	    double max_improv = - std::numeric_limits<double>::infinity();
	    double tau = 1.0;
	    
	    for (int k = 0; k != model_size; ++ k)
	      if (x[i] > 0.0) {
		const double numer = x[k] * d[k];
		const double denom = x[k] * x[k] * H(k, k);
		
		if (denom > 0.0) {
		  const double improv = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
		  
		  if (improv > max_improv) {
		    max_improv = improv;
		    tau = std::max(0.0, std::min(1.0, numer / denom));
		    v = k;
		  }
		}
	      }
	    
	    if (v  >= 0) {
	      double update = x[v] * tau;
	      if (neq + update < 0.0)
		update = - neq;
	      
	      if (update != 0.0) {
		updated = true;
		
		neq += update;
		x[v] -= update;
		
		for (int i = 0; i != model_size; ++ i)
		  d[i] -= update * H(i, v);
	      }
	    }
	  }
	  
	  if (! upated) break;
	  
	  // compute objectives
	  objective_primal = 0.0;
	  objective_dual   = 0.0;
	  for (int i = 0; i != model_size; ++ i) {
	    objective_primal += 0.5 * x[i] * (f[i] + d[i]);
	    objective_dual   += 0.5 * x[i] * (f[i] - d[i]);
	    objective_dual   += C * d[i];
	  }
	  
	  if (objective_primal - objective_dual <= tolerance) break;
	}

	return objective_primal;
      }
    };
  };
};

#endif
