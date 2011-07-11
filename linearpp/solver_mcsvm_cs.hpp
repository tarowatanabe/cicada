// -*- mode: c++ -*-

#ifndef __LINEARPP__SOLVER_MCSVMCS__HPP__
#define __LINEARPP__SOLVER_MCSVMCS__HPP__ 1

namespace linearpp
{
  // A coordinate descent algorithm for 
  // multi-class support vector machines by Crammer and Singer
  //
  //  min_{\alpha}  0.5 \sum_m ||w_m(\alpha)||^2 + \sum_i \sum_m e^m_i alpha^m_i
  //    s.t.     \alpha^m_i <= C^m_i \forall m,i , \sum_m \alpha^m_i=0 \forall i
  // 
  //  where e^m_i = 0 if y_i  = m,
  //        e^m_i = 1 if y_i != m,
  //  C^m_i = C if m  = y_i, 
  //  C^m_i = 0 if m != y_i, 
  //  and w_m(\alpha) = \sum_i \alpha^m_i x_i 
  //
  // Given: 
  // x, y, C
  // eps is the stopping tolerance
  //
  // solution will be put in w
  //
  // See Appendix of LIBLINEAR paper, Fan et al. (2008)

  class SolverMCSVMCS
  {
  public:
    
    double solve(double* w)
    {
      
    }

    void solve_sub_problem()
    {
      
      
    }
    
    int nr_class;
    
  };
};

#endif

