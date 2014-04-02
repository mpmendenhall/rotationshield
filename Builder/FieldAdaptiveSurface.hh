/* 
 * FieldAdaptiveSurface.hh, part of the RotationShield program
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

#ifndef FIELDADAPTIVESURFACE_HH
#define FIELDADAPTIVESURFACE_HH 1

#include "DVFunc.hh"
#include "SurfaceProfiles.hh"
#include "FieldEstimator2D.hh"
#include "BicubicGrid.hh"


/// Profile for optimizing cylinder to field
class FieldAdaptiveSurface: public DVFunc1<2,double> {
public:
	/// constructor
	FieldAdaptiveSurface(const DVFunc1<2,double>& f);
	
	/// evaluate function
	virtual vec2 operator()(double x) const;
	
	/// derivative
	virtual vec2 deriv(double x) const;
	
	/// optimize for field configuration
	void optimizeSpacing(const FieldEstimator2D& fes, double pfixed = 0.5, bool useDeriv = true);
	/// set constant spacing
	void setConstantSpacing();
	
	/// print out test points
	void symmetry_test() const;
	
protected:
	/// derivative of l distortion parameter
	double l_dist_deriv(double l) const;
	
	/// wrap a number into [0,1) for periodic function, starting half-way on outside
	double wrap(double x) const { double i; return x>=0 ? modf(x,&i) : 1+modf(x,&i); }
	
	CubicGrid l_remap;	//< distortion function, as interpolator
	
	const DVFunc1<2,double>& F;	//< reference function being distorted
};

#endif
