/* 
 * Interpolator.hh, part of the RotationShield program
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

#ifndef INTERPOLATOR_HH
/// make sure this file is loaded only once
#define INTERPOLATOR_HH

#include <vector>
#include <cassert>
#include <stdio.h>
#include <cmath>

/// boundary conditions for interpolation
enum BoundaryCondition {
	BC_INFINITE,		///< extend endpoint values out to infinity
	BC_CYCLIC,			///< cyclic boundary conditions
	BC_DERIVCLAMP_ZERO	///< BC for clamping derivative to 0 at 0 for bicubic interpolation
};

/// infinite, bi-directional sequence of points
class DataSequence {
public:
	/// constructor
	DataSequence(BoundaryCondition b = BC_CYCLIC): npts(0), bc(b) {}
	/// destructor
	virtual ~DataSequence() {}
	/// data retreival
	virtual double valueAt(int i, void* xopts) const = 0;
	/// get number of internal points
	virtual int getNpts() { return npts; }
	/// display state
	virtual void display() const { printf("<Generic DataSequence>\n"); }
protected:
	/// coerce data point into actual data range
	virtual unsigned int coerce(int i) const {
		if(bc == BC_CYCLIC) {
			if(i>=0)
				return i%npts;
			return ((i%npts)+npts)%npts;
		} else if(bc == BC_INFINITE) {
			if(i<=0) return 0;
			if(i>=npts) return npts-1;
			return (unsigned int)i;
		} else if (bc == BC_DERIVCLAMP_ZERO) {
			assert(npts>=2);
			if(i<0) return 1;
			if(i>=npts) return npts-1;
			return (unsigned int)i;
		}
		return 0;
	}
	int npts;				///< number of internal gridpoints
	BoundaryCondition bc;	///< boundary conditions for generating sequence
};

/// data sequence based on internal array of doubles
class DoubleSequence: public DataSequence {
public:
	DoubleSequence(BoundaryCondition b = BC_CYCLIC): DataSequence(b) {}
	/// value at given grid location
	virtual double valueAt(int i, void*) const { return pts[coerce(i)]; }
	/// add data point
	void addPoint(double p) { pts.push_back(p); npts++; }
	/// display state
	virtual void display() const {
		printf("<DoubleSequence> "); 
		for(unsigned int i = 0; i<pts.size(); i++)
			printf("%f, ",pts[i]);
		printf("(%i)\n",npts);
	}
protected:	
	std::vector<double> pts;	///< internal list of points
};

/// generic interpolator for intermediate points in a sequence
class Interpolator {
public:
	/// constructor, with input scale factor s, offset o (where interpolator '0' point is placed on input axis)
	Interpolator(DataSequence* L, double s = 1.0, double o = 0.0): scale(s), offset(o), myData(L) {}
	/// destructor
	virtual ~Interpolator() {}
	/// subclass interpolation method here; nearest-neighbor interpolation used by default
	virtual double eval(const double* x) const {
		double y;
		int i = locate(*x,&y);
		if(y<=0.5)
			return myData->valueAt(i,(void*)(x+1));
		return myData->valueAt(i+1,(void*)(x+1));
	}
	/// return interpolation basis function of distance from point
	virtual double basisFunc(double x) const { return fabs(x)<0.5; }
	/// make a new interpolator (interchangeable function pointer)
	static Interpolator* newInterpolator(DataSequence* L) { return new Interpolator(L); }
	/// set half-space offset for symmetric distribution in o0 to o0+s
	void setSymmetricOffset(double o0 = 0) { offset = o0 + 0.5*scale/myData->getNpts(); }
	
	double scale;			///< internal length scale
	double offset;			///< zero point coordinate offset
	
protected:
	/// locate position in data coordinates, remainder in [0,1)
	virtual int locate(double x, double* remainder = NULL) const {
		double l = (x-offset)*myData->getNpts()/scale;
		if(remainder)
			*remainder = l-floor(l);
		return int(floor(l));
	}
	DataSequence* myData;	///< sequence to be interpolated
};

/// linear interpolator
class LinTerpolator: public Interpolator {
public:
	/// constructor, with scale factor s, offset o
	LinTerpolator(DataSequence* L, double s = 1.0, double o = 0.0): Interpolator(L,s,o) {}
	/// evaluation by linear interpolation
	virtual double eval(const double* x) const {
		double y;
		int i = locate(*x,&y);
		return myData->valueAt(i,(void*)(x+1))*(1-y)+myData->valueAt(i+1,(void*)(x+1))*y;
	}
	/// return interpolation basis function of distance from point
	virtual double basisFunc(double x) const {
		x = fabs(x);
		if(x<1) return 1-x;
		return x;
	}
	/// make a new interpolator (interchangeable function pointer)
	static Interpolator* newLinTerpolator(DataSequence* L) { return new LinTerpolator(L); }
};



/// cubic interpolator
class CubiTerpolator: public Interpolator {
public:
	/// constructor, with scale factor s, offset o, 'sharpening' a
	CubiTerpolator(DataSequence* L, double s = 1.0, double o = 0.0, double a = -0.5): Interpolator(L,s,o), A(a) { }
	/// evaluation by linear interpolation
	virtual double eval(const double* x) const {
		double y;
		int i = locate(*x,&y);
		double p0 = myData->valueAt(i-1,(void*)(x+1));
		double p1 = myData->valueAt(i,(void*)(x+1));
		double p2 = myData->valueAt(i+1,(void*)(x+1));
		double p3 = myData->valueAt(i+2,(void*)(x+1));
		return ( A*p0*(1-y)*(1-y)*y
				+p1*(1-y)*(1-y*((2+A)*y-1))
				-p2*y*(A*(1-y)*(1-y)+y*(2*y-3))
				+A*p3*(1-y)*y*y );

	}
	/// return interpolation basis function of distance from point
	virtual double basisFunc(double x) const {
		x = fabs(x);
		if(x<1) return (1-x)*(1-x*((2+A)*x-1));
		else if(x<2) return A*(2-x)*(2-x)*(x-1);
		else return 0;
	}
	/// make a new interpolator (interchangeable function pointer)
	static Interpolator* newCubiTerpolator(DataSequence* L) { return new CubiTerpolator(L); }
	
protected:
	double A;	///< "sharpening" coefficient, default = -0.5
};

/// a sequence of interpolators for multi-dimensional interpolation between interpolators
class InterpoSequence: public DataSequence {
public:
	/// constructor
	InterpoSequence(BoundaryCondition b = BC_CYCLIC): DataSequence(b) {}
	/// add data point
	void addPoint(Interpolator* i) { myInterpolators.push_back(i); npts++; }
	/// data retreival
	virtual double valueAt(int i, void* xopts) const { return myInterpolators[coerce(i)]->eval((const double*)xopts); }
protected:
	std::vector<Interpolator*> myInterpolators;	///< interpolators in sequence
};


#endif
