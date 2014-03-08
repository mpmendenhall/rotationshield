/// \file "linmin.hh" \brief Least-squares linear polynomial fits
#ifndef LINMIN_HH
/// Make sure this header is included only once
#define LINMIN_HH 1

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include <iostream>
#include "Polynomial.hh"
#include "Typedefs.hh"

/// least-squares minimize coeffs*x = rslt+resid using QR decomposition. *frees* coeffs; needs proper-sized resid, returns x.
gsl_vector* lsmin(gsl_matrix* coeffs, const gsl_vector* rslt, gsl_vector* resid);
/// Linear fit to 3-variate polynomial; return rms residual. Coords is an N*3 matrix for N values at locations x,y,z
double polynomialFit(const gsl_matrix* coords, const gsl_vector* values, Polynomial<3,mdouble>& p);

#endif
