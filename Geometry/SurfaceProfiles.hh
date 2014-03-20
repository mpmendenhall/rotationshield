#ifndef SURFACEPROFILES_HH
#define SURFACEPROFILES_HH 1

#include "SurfaceGeometry.hh"

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
	Arc2D(double rr, vec2 v0 = vec2(), double th0 = -M_PI, double th1 = M_PI): r(rr), v(v0), t0(th0), dt(th1-th0) {}
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


/// Parameter distorter for speeding up or slowing down passage through segment in PathJoiner
class PathStretcher: public DVFunc1<2,mdouble> {
public:
	/// constructor
	PathStretcher(DVFunc1<2,mdouble>* f0, mdouble a, bool dodelete = true): u(a), f(f0), delete_f(dodelete) { assert(f); }
	/// destructor
	virtual ~PathStretcher() { if(delete_f) delete f; }
	/// evaluate function
	virtual vec2 operator()(mdouble x) const { return (*f)(distort(x)); }
	/// derivative
	virtual vec2 deriv(mdouble x) const { return f->deriv(distort(x)) * d_distort(x); }

	mdouble u;				//< distortion parameter, 1 for no distortion, 0 to infty for slow/fast stretch
	
protected:
	/// distortion map function
	virtual mdouble distort(mdouble l) const { return l * (l*(2*l-3)*(u-1) + u) ; }
	/// distortion derivative function
	virtual mdouble d_distort(mdouble l) const { return 6*l*(l-1)*(u-1) + u; }

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

/// 



#endif
