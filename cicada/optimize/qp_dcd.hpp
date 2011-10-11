// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPTIMIZE__QP_DCD__HPP__
#define __CICADA__OPTIMIZE__QP_DCD__HPP__ 1

//
// Algorithm 3 of dual coordinate descent method
//
// @inproceedings{Hsieh:2008:DCD:1390156.1390208,
// author = {Hsieh, Cho-Jui and Chang, Kai-Wei and Lin, Chih-Jen and Keerthi, S. Sathiya and Sundararajan, S.},
// title = {A dual coordinate descent method for large-scale linear SVM},
// booktitle = {Proceedings of the 25th international conference on Machine learning},
// series = {ICML '08},
// year = {2008},
// isbn = {978-1-60558-205-4},
// location = {Helsinki, Finland},
// pages = {408--415},
// numpages = {8},
// url = {http://doi.acm.org/10.1145/1390156.1390208},
// doi = {http://doi.acm.org/10.1145/1390156.1390208},
// acmid = {1390208},
// publisher = {ACM},
// address = {New York, NY, USA},
//}

//
// min   x^{\top} H x + f^{\top} * x
// 0 \leq x[i] \leq C
//

// H(i, j) = h_i^{\top} \cdot h_j

// M(w, x): w = \sum_j x_j h_j
// M(w, i): returns w^{\top} \cdot h_i
// M(w, update, i): w = w + update h_i

#include <numeric>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <cicada/weight_vector.hpp>

namespace cicada
{
  namespace optimize
  {
    
    struct QPDCD
    {
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      template <typename __X, typename __F, typename __H, typename __M>
      double operator()(__X& x,
			const __F& f,
			const __H& H,
			const __M& M,
			const double C,
			const double tolerance)
      {
	typedef cicada::WeightVector<double> weights_type;
	typedef std::vector<double, std::allocator<double> > q_type;
	typedef std::vector<size_type, std::allocator<size_type> > active_type;
	
	if (x.size() != f.size())
	  throw std::runtime_error("size does not match");
	
	const size_type model_size = f.size();
	size_type active_size = f.size();
	
	// initialize w by \sum x[i] M[i]
	weights_type w;
	M(w, x);
	
	q_type      QD(model_size, 0.0); // Q + D (we will use H notation, though... for L_1 SVM, D_ii = 0)
	active_type actives(model_size);
	for (size_type i = 0; i != model_size; ++ i) {
	  actives[i] = i;
	  QD[i] = H(i, i);
	}

	double PGmax_old =   std::numeric_limits<double>::infinity();
	double PGmin_old = - std::numeric_limits<double>::infinity();
	
	// max 1000 iterations...
	for (int iter = 0; iter != 1000; ++ iter) {
	  double PGmax_new = - std::numeric_limits<double>::infinity();
	  double PGmin_new =   std::numeric_limits<double>::infinity();
	  
	  std::random_shuffle(actives.begin(), actives.begin() + active_size);
	  
	  for (size_type s = 0; s < active_size; ++ s) {
	    const size_type i = actives[s];
	    const double G = M(w, i) + f[i];
	    
	    double PG = 0.0;
	    if (x[i] == 0.0) {
	      if (G > PGmax_old) {
		-- active_size;
		std::swap(actives[s], actives[active_size]);
		-- s;
		continue;
	      } else if (G < 0.0)
		PG = G;
	    } else if (x[i] == C) {
	      if (G < PGmin_old) {
		-- active_size;
		std::swap(actives[s], actives[active_size]);
		-- s;
		continue;
	      } else if (G > 0.0)
		PG = G;
	    } else
	      PG = G;
	    
	    PGmax_new = std::max(PGmax_new, PG);
	    PGmin_new = std::min(PGmin_new, PG);
	    
	    if (std::fabs(PG) > 1e-12) {
	      const double x_old = x[i];
	      x[i] = std::min(std::max(x[i] - G / QD[i], 0.0), C);
	      M(w, x[i] - x_old, i);
	    }
	  }
	  
	  if ((PGmax_new - PGmin_new) <= tolerance) {
	    if (active_size == model_size)
	      break;
	    else {
	      // no shrinking for the next iteration
	      active_size = model_size;
	      PGmax_old =   std::numeric_limits<double>::infinity();
	      PGmin_old = - std::numeric_limits<double>::infinity();
	      continue;
	    }
	  }
	  
	  PGmax_old = (PGmax_new <= 0.0 ?   std::numeric_limits<double>::infinity() : PGmax_new);
	  PGmin_old = (PGmin_new >= 0.0 ? - std::numeric_limits<double>::infinity() : PGmin_new);
	}
	// normalize x!
	const double sum = std::accumulate(x.begin(), x.end(), 0.0);
	if (sum != 0.0)
	  std::transform(x.begin(), x.end(), x.begin(), std::bind2nd(std::multiplies<double>(), 1.0 / sum));
	
	std::vector<double, std::allocator<double> > d(f.begin(), f.end());
	for (size_t i = 0; i != model_size; ++ i)
	  if (x[i] > 0.0)
	    for (size_t j = 0; j != model_size; ++ j)
	      d[j] += H(j, i) * x[i];
	
	double objective_primal = 0.0;
	for (size_t i = 0; i != model_size; ++ i)
	  objective_primal += 0.5 * x[i] * (f[i] + d[i]);
	
	return objective_primal;
      }
    };
  };
};

#endif
