#ifndef SURFACEPROFILES_HH
#define SURFACEPROFILES_HH 1

#include "SurfaceGeometry.hh"
#include <cmath>

/// Convenience linear sweep DFunc
class Line2D: public DVFunc1<2,mdouble> {
public:
	/// constructor
	Line2D(vec2 a, vec2 b): x0(a), x1(b) {}
	/// function call
	virtual vec2 operator()(mdouble x) const { return x0*(1-x) + x1*x; }
	/// derivative
	virtual vec2 deriv(mdouble) const { return x1-x0; }
	
	vec2 x0,x1;	//< endpoints
};

/// Arc
class Arc2D: public DVFunc1<2,mdouble> {
public:
	/// constructor
	Arc2D(double rr, double th0 = -M_PI, double th1 = M_PI, vec2 v0 = vec2()): r(rr), v(v0), t0(th0), dt(th1-th0) {}
	/// function call
	virtual vec2 operator()(mdouble x) const;
	/// derivative
	virtual vec2 deriv(mdouble x) const;
	
	double r;	//< radius
	vec2 v;		//< center
	double t0;	//< start angle
	double dt;	//< angular range
};

/// Parabolic dish
class Dish: public DVFunc1<2,mdouble> {
public:
	/// constructor
	Dish(double Z0, double DZ, double R): z0(Z0), dz(DZ), r(R) {}
	/// function call
	virtual vec2 operator()(mdouble x) const;
	
	double z0;
	double dz;
	double r;
};

/// Parameter distorter for speeding/slowing arclength traversal at ends of path
class ParamStretcher: public DVFunc1<2,mdouble> {
public:
	/// constructor
	ParamStretcher(DVFunc1<2,mdouble>* f0, mdouble d0, mdouble r0, mdouble d1, mdouble r1, bool dodelete = true);
	/// destructor
	virtual ~ParamStretcher() { if(delete_f) delete f; }
	
	/// evaluate function
	virtual vec2 operator()(mdouble x) const { return (*f)(distort(x)); }
	/// derivative
	virtual vec2 deriv(mdouble x) const { return f->deriv(distort(x)) * d_distort(x); }
	
protected:
	/// distortion map function
	virtual mdouble distort(mdouble l) const { return h + k*(l + c0*exp(-l*ri0) + c1*exp(-(1-l)*ri1) ); }
	/// distortion derivative function
	virtual mdouble d_distort(mdouble l) const { return k*(1 -ri0*c0*exp(-l*ri0) + ri1*c1*exp(-(1-l)*ri1) ); }

	// internal parameters
	double ri0;
	double ri1;
	double c0;
	double c1;
	double k;
	double h;
	
	DVFunc1<2,mdouble>* f;	//< function being distorted
	bool delete_f;			//< whether to delete function f on destruction
};

/// Circular-rounded-ends circular slab
class RoundedSlab: public PathJoiner<2,mdouble> {
public:
	/// constructor
	RoundedSlab(mdouble z0, mdouble r, mdouble w);
	/// destructor
	virtual ~RoundedSlab() { clear(); }
};

/// Rounded-end finite-thickness tube
class RoundedTube: public PathJoiner<2,mdouble> {
public:
	/// constructor
	RoundedTube(vec2 x0, vec2 x1, mdouble r, mdouble endfrac = 0.33);
	/// destructor
	virtual ~RoundedTube() { clear(); }
	
	/// evaluate function
	virtual vec2 operator()(mdouble x) const { return PathJoiner<2,mdouble>::operator()(wrap(x)); }
	/// derivative
	virtual vec2 deriv(mdouble x) const { return PathJoiner<2,mdouble>::deriv(wrap(x)); }

	/// wrap a number into [0,1] for periodic function
	mdouble wrap(mdouble x) const { mdouble i; x += 0.5*seg_endpts[1]/seg_endpts.back(); return x>0 ? modf(x,&i) : 1-modf(x,&i); }
};



#endif
