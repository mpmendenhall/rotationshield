/* 
 * Tests.hh, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

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
void superball_test(unsigned int ngrid = 32);
/// endless costheta between SC caps
void mirror_test();
/// thick ferromagnetic tube
void tube_test();
/// calculate trapped-flux state of superconducting loop
void flux_trap_test();

#endif
