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
//
// sum x = C
// x \geq 0
//

// H(i, j) = h_i^{\top} \cdot h_j

// M(w, x): w = \sum_j x_j h_j
// M(w, i): returns w^{\top} \cdot h_i
// M(w, update, i): w = w + update h_i

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
      
      template <typename X, typename F, typename H, typename M>
      double operator()(X& x,
			const F& f,
			const H& H,
			const M& M,
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
	  index[i] = i;
	  QD[i] = H(i, i);
	}

	double PGmax_old =   std::numeric_limits<double>::infinity();
	double PGmin_old = - std::numeric_limits<double>::infinity();
	
	// max 1000 iterations...
	for (int iter = 0; iter != 1000; ++ iter) {
	  double PGmax_new = - std::numeric_limits<double>::infinity();
	  double PGmin_new =   std::numeric_limits<double>::infinity();
	  
	  std::random_shuffle(actives.begin(), actives.begin() + active_size);
	  
	  active_type::iterator aiter_end = actives.end();
	  for (active_type::iterator aiter = actives.begin(); aiter != aiter_end; ++ aiter) {
	    const size_type i = *aiter;
	    const double G = M(weights, i) + f[i];
	    
	    double PG = 0.0;
	    if (x[i] == 0.0) {
	      if (G > PGmax_old) {
		-- active_size;
		swap(*aiter, index[active_size]);
		-- aiter;
		continue;
	      } else if (G < 0.0)
		PG = G;
	    } else if (x[i] == C) {
	      if (G < PGmin_old) {
		-- active_size;
		swap(*aiter, index[active_size]);
		-- aiter;
		continue;
	      } else if (G > 0.0)
		PG = G;
	    } else
	      PG = G;
	    
	    PGmax_new = max(PGmax_new, PG);
	    PGmin_new = min(PGmin_new, PG);
	    
	    if (std::fabs(PB) > 1e-12) {
	      const dobule x_old = x[i];
	      x[i] = std::min(std::max(x[i] - G / QD[i]), C);
	      M(weights, x[i] - x_old, i)
	    }
	  }
	  
	  if ((PGmax_new - PGmin_new) <= toleranse) {
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
      }
    };
  };
};

#endif
