/* 
 * SurfaceProfiles.hh, part of the RotationShield program
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

#ifndef SURFACEPROFILES_HH
#define SURFACEPROFILES_HH 1

#include "SurfaceGeometry.hh"
#include <cmath>

/// Convenience linear sweep DFunc
class Line2D: public DVFunc1<2,double> {
public:
	/// constructor
	Line2D(vec2 a, vec2 b): x0(a), x1(b) {}
	/// function call
	virtual vec2 operator()(double x) const { return x0*(1-x) + x1*x; }
	/// derivative
	virtual vec2 deriv(double, bool normalized = false) const { if(normalized) return (x1-x0).normalized(); return x1-x0; }
	
	vec2 x0,x1;	//< endpoints
};

/// Arc of circle
class Arc2D: public DVFunc1<2,double> {
public:
	/// constructor
	Arc2D(double rr, double th0 = -M_PI, double th1 = 0, vec2 v0 = vec2());
	/// function call
	virtual vec2 operator()(double x) const;
	/// derivative
	virtual vec2 deriv(double x, bool normalized = false) const;
	
	double r;	//< radius
	vec2 v;		//< center
	double t0;	//< start angle
	double dt;	//< angular range
};

/// Parabolic dish
class Dish: public DVFunc1<2,double> {
public:
	/// constructor
	Dish(double Z0, double DZ, double R): z0(Z0), dz(DZ), r(R) {}
	/// function call
	virtual vec2 operator()(double x) const;
	
	double z0;
	double dz;
	double r;
};

/// line with cosine wiggles
class CosineLine: public DVFunc1<2,double> {
public:
	/// constructor
	CosineLine(vec2 start, vec2 end, unsigned int ncyc, double ampl);
	/// function call
	virtual vec2 operator()(double x) const;
	
protected:
	vec2 x0,x1,o;	//< endpoints and orthogonal
	unsigned int n;	//< number of cycles
};

/// Parameter distorter for speeding/slowing arclength traversal at ends of path
class ParamStretcher: public DVFunc1<2,double> {
public:
	/// constructor
	ParamStretcher(DVFunc1<2,double>* f0, double d0, double r0, double d1, double r1, bool dodelete = true);
	/// destructor
	virtual ~ParamStretcher() { if(delete_f) delete f; }
	
	/// evaluate function
	virtual vec2 operator()(double x) const { return (*f)(distort(x)); }
	/// derivative
	virtual vec2 deriv(double x, bool normalized = false) const { return f->deriv(distort(x),normalized) * (normalized?1:d_distort(x)); }
	
protected:
	/// distortion map function
	virtual double distort(double l) const { return h + k*(l + c0*exp(-l*ri0) + c1*exp(-(1-l)*ri1) ); }
	/// distortion derivative function
	virtual double d_distort(double l) const { double dd = k*(1 -ri0*c0*exp(-l*ri0) + ri1*c1*exp(-(1-l)*ri1) );  assert(dd); return dd; }

	// internal parameters
	double ri0;
	double ri1;
	double c0;
	double c1;
	double k;
	double h;
	
	DVFunc1<2,double>* f;	//< function being distorted
	bool delete_f;			//< whether to delete function f on destruction
};

//
//==========================================================================
//


/// append sequential univariate functions, assumed to each reside on [0,1]
template<unsigned int N, typename T>
class PathJoiner: public DVFunc1<N,T> {
public:
	/// constructor
	PathJoiner() {
		seg_endpts.push_back(0);
	}
	
	/// destructor
	virtual ~PathJoiner() {}
	
	/// evaluate function
	virtual Vec<N,T> operator()(T x) const {
		if(this->period) x = wrap(x);
		unsigned int i = locate(x);
		return seg_offsets[i] + (*mysegs[i])(x);
	}
	
	/// derivative
	virtual Vec<N,T> deriv(T x, bool normalized = false) const {
		if(this->period) x = wrap(x);
		unsigned int i = locate(x);
		return mysegs[i]->deriv(x,normalized) * (normalized?1:seg_endpts.back()/(seg_endpts[i+1]-seg_endpts[i]));
	}
	
	/// append a segment, spaced to maintain continuous pathlength derivative
	virtual void append( DVFunc1<N,T>* f, bool smooth_pathlength = true) {
		assert(f);
		if(!mysegs.size()) {
			seg_endpts.push_back(1);
		} else {
			unsigned int i = mysegs.size();
			T mul = smooth_pathlength ? f->deriv(0).mag() / mysegs.back()->deriv(1).mag() * (seg_endpts[i]-seg_endpts[i-1]) : 1;
			seg_endpts.push_back(seg_endpts.back() + mul);
		}
		if(mysegs.size())
			seg_offsets.push_back(seg_offsets.back() + (*mysegs.back())(1.) - (*f)(0));
		else
			seg_offsets.push_back(Vec<N,T>());
		mysegs.push_back(f);
		if( ((*mysegs.back())(1)+seg_offsets.back()-(*mysegs[0])(0)).mag2() < 1e-6 ) this->period = 1;
		else this->period = 0;
	}
	
	/// display contents
	void display() const {
		std::cout << "PathJoiner in " << mysegs.size() << " segments:\n";
		for(unsigned int i=0; i<mysegs.size(); i++)
			std::cout << "\t" << seg_endpts[i] << " to " << seg_endpts[i+1] << "\tstarting at " << seg_offsets[i] << std::endl;
	}
	
protected:
	
	/// locate sub-function and position therein
	unsigned int locate(T& l) const {
		assert(mysegs.size());
		typename std::vector<T>::const_iterator it = std::upper_bound(seg_endpts.begin(), seg_endpts.end(), l*seg_endpts.back());
		unsigned int i = it-seg_endpts.begin();
		if(!i) i += 1;
		if(it==seg_endpts.end()) i -= 1;
		l = (l*seg_endpts.back()-seg_endpts[i-1])/(seg_endpts[i]-seg_endpts[i-1]);
		return i-1;
	}
	
	/// delete all items
	void clear() {
		for(unsigned int i=0; i<mysegs.size(); i++)
			delete mysegs[i];
		mysegs.clear();
		seg_endpts.clear();
		seg_offsets.clear();
		seg_endpts.push_back(0);
	}
	
	/// wrap a number into [0,1) for periodic function, starting half-way on outside
	T wrap(T x) const { T i; x += 0.5*seg_endpts[1]/seg_endpts.back(); return x>=0 ? modf(x,&i) : 1+modf(x,&i); }
	
	std::vector< DVFunc1<N,T>* > mysegs;	//< subsections
	std::vector<T> seg_endpts;				//< map from global l to sub-range endpoints
	std::vector< Vec<N,T> > seg_offsets;	// endpoint offsets for each segment
};

/// Circular-rounded-ends circular slab
class RoundedSlab: public PathJoiner<2,double> {
public:
	/// constructor
	RoundedSlab(double z0, double r, double w, double endfrac = 0.33);
	/// destructor
	virtual ~RoundedSlab() { clear(); }
};

/// Slab with CosineLine wiggly sides
class WiggleSlab: public PathJoiner<2,double> {
public:
	/// constructor
	WiggleSlab(double z0, double r, double w, unsigned int n, double a, double endfrac = 0.33);
	/// destructor
	virtual ~WiggleSlab() { clear(); }
};


/// Rounded-end finite-thickness tube
class RoundedTube: public PathJoiner<2,double> {
public:
	/// constructor
	RoundedTube(vec2 x0, vec2 x1, double r, double endfrac = 0.33);
	/// destructor
	virtual ~RoundedTube() { clear(); }
};



#endif
