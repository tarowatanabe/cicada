/* The operation of the code is mostly controlled by the parameters
   in the cg_parameter structure.  In the following example,
   the parameter QuadStep is set to FALSE.  When QuadStep is TRUE,
   the trial step in each iteration is computed as the minimizer
   of a quadratic interpolant along the search direction. In
   performing the quadstep, we hope to find a suitable line search
   point right away, completely by-passing the secant iteration.
   However, as the iterates approach a minimizer, the numerical accuracy of
   the minimizer of the quadratic interpolant becomes worse. When the relative
   change in the function values for two consecutive iterations reach
   QuadCutOff, then the code completely turns off the quadstep. The user
   can turn off the quadstep by setting QuadStep FALSE. By leaving
   QuadStep TRUE, but increasing QuadCutOff (default 1.d-12), the code turns
   off the QuadStep sooner. Below, we run the code twice, first with
   the QuadStep turned off, then with the QuadStep turned on. Notice that
   the performance improves with the QuadStep is on.  This behavior is typical.

   Termination status: 0
   Convergence tolerance for gradient satisfied
   maximum norm for gradient:  4.823094e-09
   function value:            -6.530787e+02

   cg  iterations:                  32
   function evaluations:            36
   gradient evaluations:            68
   ===================================

   Termination status: 0
   Convergence tolerance for gradient satisfied
   maximum norm for gradient:  6.269565e-09
   function value:            -6.530787e+02

   cg  iterations:                  32
   function evaluations:            54
   gradient evaluations:            46 */

#include <math.h>
#include "cg_user.h"

double myvalue
(
    double   *x,
    INT       n
) ;

void mygrad
(
    double    *g,
    double    *x,
    INT        n
) ;

double myvalgrad
(
    double    *g,
    double    *x,
    INT        n
) ;

int main (void)
{
    double *x ;
    INT i, n ;
    cg_parameter Parm ;

    /* allocate space for solution */
    n = 100 ;
    x = (double *) malloc (n*sizeof (double)) ;

    /* set starting guess */
    for (i = 0; i < n; i++) x [i] = 1. ;

    cg_default (&Parm) ;    /* set default parameter values */
    Parm.QuadStep = FALSE ; /* change QuadStep to FALSE */

    /* run the code */
    cg_descent(x, n, NULL, &Parm, 1.e-8, myvalue, mygrad, myvalgrad, NULL) ;

    /* set starting guess */
    for (i = 0; i < n; i++) x [i] = 1. ;
    Parm.QuadStep = TRUE ; /* change QuadStep to TRUE */

    /* run the code */
    cg_descent(x, n, NULL, &Parm, 1.e-8, myvalue, mygrad, myvalgrad, NULL) ;

    free (x) ; /* free workspace */
}

double myvalue
(
    double   *x,
    INT       n
)
{
    double f, t ;
    INT i ;
    f = 0. ;
    for (i = 0; i < n; i++)
    {
        t = i+1 ;
        t = sqrt (t) ;
        f += exp (x [i]) - t*x [i] ;
    }
    return (f) ;
}

void mygrad
(
    double    *g,
    double    *x,
    INT        n
)
{
    double t ;
    INT i ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        g [i] = exp (x [i]) -  t ;
    }
    return ;
}

double myvalgrad
(
    double    *g,
    double    *x,
    INT        n
)
{
    double ex, f, t ;
    INT i ;
    f = (double) 0 ;
    for (i = 0; i < n; i++)
    {
        t = i + 1 ;
        t = sqrt (t) ;
        ex = exp (x [i]) ;
        f += ex - t*x [i] ;
        g [i] = ex -  t ;
    }
    return (f) ;
}
