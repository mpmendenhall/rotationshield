/* 
 * InterpolatingRS.hh, part of the RotationShield program
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

/// \file InterpolatingRS.hh \brief Virtual base class for interpolated collections of interacting degrees of freedom

#ifndef INTERPOLATINGRS_HH
/// Make sure this header is only loaded once
#define INTERPOLATINGRS_HH

#include "ReactiveSet.hh"
#include "InterpolationHelper.hh"
#include "BicubicGrid.hh"
#include <vector>

/// ReactiveSet with utilities for interpolating between DF
class InterpolatingRS: public ReactiveSet {
public:
	/// constructor
	InterpolatingRS(unsigned int nph): ReactiveSet(nph) {}
	
	//===================================== ReactiveSet subclass
	/// total number of degrees of freedom
	virtual unsigned int nDF() const { return InterplDF.n_pts(); }
	//=====================================
	
protected:

	//===================================== ReactiveSet subclass
	/// additional routines for setting a DF value
	virtual void _setDF(unsigned int DF, double v) { InterplDF[DF] = v; }
	/// additional routines for setting entire state vector
	virtual void _setDFv(const mvec& v) { assert(v.size() == nDF()); InterplDF.setData(&v[0]); }
	//=====================================
	
	InterpolationHelper InterplDF;	///< degrees of freedom in InterpolatingHelper grid
};

/// Simplified case for 2D interpolation
class InterpolatingRS2D: public ReactiveSet {
public:
	/// constructor
	InterpolatingRS2D(unsigned int nph): ReactiveSet(nph), nZ(0), nDFi(0) { }
	/// destructor
	virtual ~InterpolatingRS2D() { clear_data(); }
	
	//===================================== ReactiveSet subclass
	/// total number of degrees of freedom
	virtual unsigned int nDF() const { return nZ*nPhi*nDFi; }
	//=====================================
	
	/// evaluate interpolators at point
	mvec interpl_DF(vec2 l) const;
	
	/// display data grids
	virtual void printData() const;
	
protected:
	
	//===================================== ReactiveSet subclass
	/// additional routines for setting a DF value
	virtual void _setDF(unsigned int DF, double v);
	//=====================================
	
	/// determine "address" for DF in interpolating grids
	void DF_address(unsigned int DF, unsigned int& p, unsigned int& z, unsigned int& d) const;
	/// delete previous grids
	virtual void clear_data();
	/// set up data grid
	virtual void make_grids(unsigned int nz, unsigned int ndf);
	
	std::vector<BicubicGrid*> G;	///< degrees of freedom stored in interpolator grid
	unsigned int nZ;				///< grid size in z direction
	unsigned int nDFi;				///< number of DF per element
};


#endif
