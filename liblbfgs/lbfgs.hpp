// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __LIBLBFGS__LBFGS__HPP__
#define __LIBLBFGS__LBFGS__HPP__ 1

#include <climits>
#include <vector>

#include <liblbfgs/lbfgs.h>
#include <liblbfgs/lbfgs_error.hpp>

namespace liblbfgs
{
  
  template <typename Function>
  class LBFGS
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    LBFGS(const Function& function,
	  const size_type max_iterations=0,
	  const lbfgsfloatval_t l1 = 0.0,
	  const size_type l1_start = 0)
      : function_(function)
    {
      lbfgs_parameter_init(&param_);
      
      param_.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
      
      if (l1 > 0.0) {
	param_.orthantwise_c = l1;
	param_.orthantwise_start = l1_start;
      } else
	param_.orthantwise_c = 0.0;
      
      param_.max_iterations = max_iterations;
    }
    
    lbfgsfloatval_t operator()(const size_type n, lbfgsfloatval_t* x)
    {
      objective_     = std::numeric_limits<lbfgsfloatval_t>::infinity();
      objective_opt_ = std::numeric_limits<lbfgsfloatval_t>::infinity();
      
      x_opt_.clear();
      x_opt_.reserve(n);
      x_opt_.resize(n);
      x_opt_.insert(x_opt_.end(), x, x + n);
      
      lbfgs(n, x, &objective_, _evaluate, 0, this, &param_);
      
      // copy the result back to x...
      std::copy(x_opt_.begin(), x_opt_.end(), x);
      
      return objective_opt_;
    }
    
  private:
    
    static lbfgsfloatval_t _evaluate(void *instance,
				     const lbfgsfloatval_t *x,
				     lbfgsfloatval_t *g,
				     const int n,
				     const lbfgsfloatval_t step)
    {
      LBFGS<Function>* obj = reinterpret_cast<LBFGS<Function>*>(instance);
      
      const lbfgsfloatval_t result = obj->function_(n, x, g);
      
      lbfgsfloatval_t result_orthantwise = result;
      if (obj->param_.orthantwise_c > 0.0) {
	lbfgsfloatval_t l1 = 0.0;
	for (difference_type i = obj->param_.orthantwise_start; i != n; ++ i)
	  l1 += std::fabs(x[i]);
	
	result_orthantwise += obj->param_.orthantwise_c * l1;
      }
      
      // keep the optimum.... this is simply a workaround...
      if (result_orthantwise < obj->objective_opt_) {
	obj->objective_opt_ = result_orthantwise;
	std::copy(x, x + n, obj->x_opt_.begin());
      }
      
      return result;
    }
    
  private:
    lbfgs_parameter_t param_;
    const Function& function_;
    
    lbfgsfloatval_t objective_;
    lbfgsfloatval_t objective_opt_;
    std::vector<lbfgsfloatval_t, std::allocator<lbfgsfloatval_t> > x_opt_;
  };
  
};

#endif
