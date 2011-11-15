#include <vector>
#include <iostream>

extern "C" {
  extern void lbfgs_fortran(int* n, int* m, double* x, double* f, double* g, 
			    int* diagco, double* diag, int* iprint, double* eps, 
			    double* xtol, double* w, int* iflag);

  extern void lbfgsb_fortran(int *n, int *m, double *x, double *l, double *u, int *nbd, 
			     double *f, double *g, double *factr, double *pgtol, double *wa, int *iwa, 
			     char *task, int *iprint, char *csave, int *lsave, int *isave, double *dsave, short task_len, short csave_len); 
}

using namespace std;

class LBFGS {

private:
  int n;
  int m;
  int iflag;
  double *l;
  double *u;
  int    *nbd;
  int type; 
  char   task[64];
  char   csave[64];
  int    lsave[4];
  int    isave[64];
  double dsave[32];

  std::vector <double> diag;
  std::vector <double> w;
  std::vector <int>    iw;

public:

  // these valuables are fixed.
  LBFGS(): n(0), m(5), iflag(0), l(0), u(0), nbd(0), type(-1) {};  
  ~LBFGS() {};

  const char *what () { return task; }

  int init (int _n, int _m, double *_l = 0, double *_u = 0, int *_nbd = 0)
  {
    n = _n;
    m = _m;

    if (_l && _u && _nbd) {
      l = _l;
      u = _u;
      nbd = _nbd;
      w.resize ((2 * m + 4) * n + 12 * m * m + 12 * m);
      iw.resize (3 * n);
      iflag = 0;
      type = 1; // LBFGS-B
    } else {
      iflag = 0;
      w.resize (n * (2 * m + 1) + 2 * m);
      diag.resize (n);
      type = 0; // LBFGS
    }

    return 0;
  }
  
  int optimize (double *x, double *f, double *g)   // vector, object, gradient
  {
    switch (type) {

    case 0: {
      int iprint[] = {-1, 0};
      double eta  = 1e-5;
      double xtol = 1e-5;
      int diagco = 0;
      lbfgs_fortran(&n, &m, x, f, g, &diagco, &diag[0], iprint, &eta, &xtol, &w[0], &iflag);
      if (iflag < 0) {
	strcpy (task, "LBFGS: routine stops with unexpected error.");
	return -1;
      } 
      if (iflag == 0) {
	strcpy (task, "LBFGS: converge");
	return 0; // terminate
      } 
      return 1; // evaluate next f and g
    }
      break;
      
    case 1: {
      int iprint = -1;
      double factr = 1.0; // may be enough
      double pgtol = 0;
      if (iflag == 0) strcpy (task, "START"); // init
      iflag = 1;

      while (true) {
	lbfgsb_fortran(&n, &m, x, l, u, nbd, f, g, &factr, &pgtol, &w[0], &iw[0], 
		task, &iprint, csave, lsave, isave, dsave, 60, 60);
	if (strncmp (task, "FG", 2) == 0) return 1;         // evaluate next f and g
	else if (strncmp (task, "NEW_X", 4) == 0) continue; // next iteration
	else if (strncmp (task, "CONV",  4) == 0) return 0; // terminate
	else return -1; // error
      }

      return 1;
    }
      break;

    default:
      strcpy (task, "LBFGS::run() call LBFGS::init() first");
      return -1;
    }

    return -1;
  }
};
