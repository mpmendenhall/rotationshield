/// \file "linmin.hh" \brief Various utilities for linear fits using CLHEP matrices
#ifndef LINMIN_HH
/// Make sure this header is included only once
#define LINMIN_HH 1

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/DiagMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include <iostream>
#include "Polynomial.hh"
#include "MiscUtils.hh"

/// Least-squares minimization
CLHEP::HepVector lsmin(CLHEP::HepMatrix coeffs, CLHEP::HepVector rslt);
/// Sum of squares of residuals from a fit
double ssresids(CLHEP::HepMatrix coeffs, CLHEP::HepVector x, CLHEP::HepVector y);
/// Multivariate quadratic least-squares fit
CLHEP::HepVector polyQuadraticFit(CLHEP::HepMatrix coords, CLHEP::HepVector values, double& resids);
/// Linear fit to 3-variate polynomial
double polynomialFit(CLHEP::HepMatrix coords, CLHEP::HepVector values, Polynomial<3,mdouble>& p);
/// Example of using polyQuadraticFit()
void test_fitting();

#endif
