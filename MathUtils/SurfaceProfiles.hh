#ifndef SURFACEPROFILES_HH
#define SURFACEPROFILES_HH 1

#include "SurfaceGeometry.hh"

/// Convenience linear sweep DFunc
class Line2D: public DVFunc1<2,mdouble> {
public:
	Line2D(vec2 a, vec2 b): x0(a), x1(b) {}
	vec2 x0,x1;
	virtual vec2 operator()(mdouble x) const { return x0*(1-x) + x1*x; }
	virtual vec2 deriv(mdouble) const { return x1-x0; }
};

/// Half-circle
class Ball: public DVFunc1<2,mdouble> {
public:
	Ball(double rr, double a = 1, double zz = 0): r(rr), z(zz), ar(a) {}
	virtual vec2 operator()(mdouble x) const { return vec2(z-r*cos(M_PI*x)*ar, r*sin(M_PI*x)); }
	double r;
	double z;
	double ar;
};

/// Parabolic dish
class Dish: public DVFunc1<2,mdouble> {
public:
	Dish(double Z0, double DZ, double R): z0(Z0), dz(DZ), r(R) {}
	virtual vec2 operator()(mdouble x) const {
		if(r>0) return vec2(z0-dz*((1-x)*(1-x)-1), r*(1-x));
		return vec2(z0-dz*(x*x-1), r*x);
	}
	
	double z0;
	double dz;
	double r;
};

class RoundedSlab: public PathJoiner<2,mdouble> {
public:
	/// constructor
	RoundedSlab(mdouble z0, mdouble r, mdouble w): PathJoiner<2,mdouble>() {
		Line2D* L1 = new Line2D(vec2(z0-0.5*w, 0), vec2(z0-0.5*w, r-0.5*w));
		Ball* B = new Ball(0.5*w);
		Line2D* L2 = new Line2D(vec2(z0+0.5*w, r-0.5*w), vec2(z0+0.5*w, 0));
		append(L1);
		append(B);
		append(L2);
	}
	
	/// destructor
	virtual ~RoundedSlab() { clear(); }
	
protected:
	mdouble z0;
	mdouble w;
	mdouble r;
};


#endif
