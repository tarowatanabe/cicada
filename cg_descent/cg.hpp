// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CG_DESCENT__CG__HPP__
#define __CG_DESCENT__CG__HPP__ 1

#include <climits>
#include <vector>

#include <cg_descent/cg_user.h>

namespace cg
{
  template <typename Function>
  class CG
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

  public:    
    CG(const Function& function,
       const size_type max_iterations=0)
      : function_(function)
    {
      cg_default(&param_);
      
      param_.maxit = max_iterations;
    }
    
    double operator()(const size_type n, double* x)
    {
      cg_descent(x, n, &stats_, &param_, 1e-5, NULL, NULL, NULL, _evaluate, NULL);
      
      return stats_.f;
    }
    
  private:
    static double _evaluate(void* instance, double* g, double* x, CG_INT n)
    {
      return reinterpret_cast<CG<Function>*>(instance)->function_(n, x, g);
    }
    
  private:
    cg_parameter param_;
    cg_stats     stats_;
    const Function& function_;
  };
};

#endif
