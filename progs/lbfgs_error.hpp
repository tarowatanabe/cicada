//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __LBFGS__ERROR_HPP__
#define __LBFGS__ERROR_HPP__ 1

#include <lbfgs.h>

static inline
const char* lbfgs_error(const int error)
{
  switch (result) {
  case LBFGS_SUCCESS:
    return "L-BFGS reaches convergence.";
  case LBFGS_ALREADY_MINIMIZED:
    return "The initial variables already minimize the objective function.";
  case LBFGSERR_UNKNOWNERROR:
    return "Unknown error.";
  case LBFGSERR_LOGICERROR:
    return "Logic error.";
  case LBFGSERR_OUTOFMEMORY:
    return "Insufficient memory.";
  case LBFGSERR_CANCELED:
    return "The minimization process has been canceled.";
  case LBFGSERR_INVALID_N:
    return "Invalid number of variables specified.";
  case LBFGSERR_INVALID_N_SSE:
    return "Invalid number of variables (for SSE) specified.";
  case LBFGSERR_INVALID_X_SSE:
    return "The array x must be aligned to 16 (for SSE).";
  case LBFGSERR_INVALID_EPSILON:
    return "Invalid parameter lbfgs_parameter_t::epsilon specified.";
  case LBFGSERR_INVALID_TESTPERIOD:
    return "Invalid parameter lbfgs_parameter_t::past specified.";
  case LBFGSERR_INVALID_DELTA:
    return "Invalid parameter lbfgs_parameter_t::delta specified.";
  case LBFGSERR_INVALID_LINESEARCH:
    return "Invalid parameter lbfgs_parameter_t::linesearch specified.";
  case LBFGSERR_INVALID_MINSTEP:
    return "Invalid parameter lbfgs_parameter_t::min_step specified.";
  case LBFGSERR_INVALID_MAXSTEP:
    return "Invalid parameter lbfgs_parameter_t::max_step specified.";
  case LBFGSERR_INVALID_FTOL:
    return "Invalid parameter lbfgs_parameter_t::ftol specified.";
  case LBFGSERR_INVALID_WOLFE:
    return "Invalid parameter lbfgs_parameter_t::wolfe specified.";
  case LBFGSERR_INVALID_GTOL:
    return "Invalid parameter lbfgs_parameter_t::gtol specified.";
  case LBFGSERR_INVALID_XTOL:
    return "Invalid parameter lbfgs_parameter_t::xtol specified.";
  case LBFGSERR_INVALID_MAXLINESEARCH:
    return "Invalid parameter lbfgs_parameter_t::max_linesearch specified.";
  case LBFGSERR_INVALID_ORTHANTWISE:
    return "Invalid parameter lbfgs_parameter_t::orthantwise_c specified.";
  case LBFGSERR_INVALID_ORTHANTWISE_START:
    return "Invalid parameter lbfgs_parameter_t::orthantwise_start specified.";
  case LBFGSERR_INVALID_ORTHANTWISE_END:
    return "Invalid parameter lbfgs_parameter_t::orthantwise_end specified.";
  case LBFGSERR_OUTOFINTERVAL:
    return "The line-search step went out of the interval of uncertainty.";
  case LBFGSERR_INCORRECT_TMINMAX:
    return "A logic error occurred; alternatively, the interval of uncertainty became too small.";
  case LBFGSERR_ROUNDING_ERROR:
    return "A rounding error occurred; alternatively, no line-search step satisfies the sufficient decrease and curvature conditions.";
  case LBFGSERR_MINIMUMSTEP:
    return "The line-search step became smaller than lbfgs_parameter_t::min_step.";
  case LBFGSERR_MAXIMUMSTEP:
    return "The line-search step became larger than lbfgs_parameter_t::max_step.";
  case LBFGSERR_MAXIMUMLINESEARCH:
    return "The line-search routine reaches the maximum number of evaluations.";
  case LBFGSERR_MAXIMUMITERATION:
    return "The algorithm routine reaches the maximum number of iterations.";
  case LBFGSERR_WIDTHTOOSMALL:
    return "Relative width of the interval of uncertainty is at most lbfgs_parameter_t::xtol.";
  case LBFGSERR_INVALIDPARAMETERS:
    return "A logic error (negative line-search step) occurred.";
  case LBFGSERR_INCREASEGRADIENT:
    return "The current search direction increases the objective function value.";
  default: return "unknown error";
  }
}

#endif
