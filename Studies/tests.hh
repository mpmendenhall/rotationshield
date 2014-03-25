/// \file "tests.hh" \brief Tests to make sure everything is working properly
#ifndef TESTS_HH
/// Make sure this header is only loaded once
#define TESTS_HH 1

#include "Typedefs.hh"

/// Print a comparison between expected and calculated values for a data point
bool compareResults(double a, double b, const char* label = 0x0);

/// Test of line source fields and symmetric solver using crudely-gridded "old-style" shield
bool reference_simpleshield();

/// tests for integrator routines
bool integrator_tests();
/// test cacheing solver to file
bool reference_simpleshield_cached();

/// superconducting ball in uniform field
void superball_test();
/// endless costheta between SC caps
void mirror_test();
/// thick ferromagnetic tube
void tube_test();
/// calculate trapped-flux state of superconducting loop
void flux_trap_test();

#endif
