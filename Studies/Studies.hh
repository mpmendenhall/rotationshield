/// \file Studies.hh \brief Various shield/coil geometry studies

#ifndef STUDIES_HH
/// Makes sure to only load this file once
#define STUDIES_HH 1

#include "CosThetaBuilder.hh"
#include "SymmetricSolver.hh"
#include "PlaneSource.hh"
#include "ShieldBuilder.hh"
#include "analysis.hh"

/// Full-scale coil
/// N=30, radius=0.648, length=4.292,
/// distortion a = -0.00295
/// shield radius nominal 6.9cm outside coil
/// shield 40cm longer than coil
void fullScaleCoil(std::ostream& gridf, std::ostream& fitf);

/// Bare Full-scale coil
/// N=30, radius=0.648, length=4.292,
/// distortion a = -0.00295
void bareFullScaleCoil(std::ostream& gridf, std::ostream& fitf);

/// Short ("vertical") coil with re-oriented sample cell
/// N=30, length=2.m, radius variable
/// shield radius nominal 10cm outside coil; shield 40cm longer than coil
/// cell dimensions y=40cm, x=7.5cm, z=10cm, inner edge 5cm from center in x (field) direction
/// Want:
///   B0 = 30mG; <dBi/di> < 10^-7 G/cm
///   sqrt(<(dBx/dz)^2>) (z=long direction, x=along polarization)
///   cells: | 7.5 | 10 | 7.5 |    --> x
void shortCoil(std::ostream& gridf, std::ostream& fitf, double rd = 1.0);

#endif
