/// \file "tests.hh" \brief Tests to make sure everything is working properly
#ifndef TESTS_HH
/// Make sure this header is only loaded once
#define TESTS_HH 1

#include "SymmetricSolver.hh"
#include "GenericSolver.hh"
#include "FieldSource.hh"
#include "PlaneSource.hh"
#include "Boxel.hh"
#include "ShieldBuilder.hh"
#include "analysis.hh"
#include "CosThetaBuilder.hh"

/// Print a comparison between expected and calculated values for a data point
bool compareResults(mdouble a, mdouble b, const char* label = 0x0);

/// Basic test that everything is working well
bool reference_sanity_check();

/// Faster test with a very crudely gridded shield
bool reference_simpleshield();

#endif
