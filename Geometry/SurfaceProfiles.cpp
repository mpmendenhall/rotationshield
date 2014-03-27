#include "SurfaceProfiles.hh"
#include <cmath>

Arc2D::Arc2D(double rr, double th0, double th1, vec2 v0): r(rr), v(v0), t0(th0), dt(th1-th0) {
	assert(dt);
	period = fabs(fmod(dt,2*M_PI)) < 1e-6;
}

vec2 Arc2D::operator()(double x) const {
	double th = t0+dt*x;
	return v + vec2(r*cos(th), r*sin(th));
}
vec2 Arc2D::deriv(double x) const {
	double th = t0+dt*x;
	return vec2(-r*sin(th)*dt, r*cos(th)*dt);
}

//----

vec2 Dish::operator()(double x) const {
	if(r>0) return vec2(z0-dz*((1-x)*(1-x)-1), r*(1-x));
	return vec2(z0-dz*(x*x-1), r*x);
}

//----

ParamStretcher::ParamStretcher(DVFunc1<2,double>* f0, double d0, double r0, double d1, double r1, bool dodelete):
ri0(1./fabs(r0)), ri1(1./fabs(r1)), c0(-(d0-1)*fabs(r0)), c1((d1-1)*fabs(r1)), k(1), h(0), f(f0), delete_f(dodelete) {
	k = 1./(distort(1)-distort(0));
	h = -distort(0);
}

//----

RoundedSlab::RoundedSlab(double z0, double r, double w, double endfrac): PathJoiner<2,double>() {

	w = copysign(w,r);
	
	assert(fabs(r) > fabs(0.5*w));
	
	double sidelen = 2*fabs(r-0.5*w);
	double endlen = fabs(0.5*w)*M_PI;
	double s = (1-endfrac)/endfrac * endlen/sidelen;

	append(new ParamStretcher(new Line2D(vec2(z0+0.5*w, 0), vec2(z0+0.5*w, fabs(r-0.5*w))),
							  1, 1, s, r));
				
	append(new Arc2D(0.5*w, 0, r<0 ? -M_PI:M_PI));
	
	append(new ParamStretcher(new Line2D(vec2(z0-0.5*w, fabs(r-0.5*w)), vec2(z0-0.5*w, 0)),
							 	s, r, 1, 1));
}

//----

RoundedTube::RoundedTube(vec2 x0, vec2 x1, double r, double endfrac): PathJoiner<2,double>() {

	assert((x1-x0).mag() > 2*fabs(r));
	assert(0 < endfrac && endfrac < 1);
	
	vec2 dl = (x1-x0).normalized()*fabs(r);
	vec2 dn = rhOrtho(dl);
	double th = atan2(dn[1],dn[0]);
	
	x1 -= dl;
	x0 += dl;
	double sidelen = (x1-x0).mag();
	double endlen = fabs(r)*M_PI;
	double s = (1-endfrac)/endfrac * endlen/sidelen;
	
	if(r>0) {
		append(new ParamStretcher(new Line2D(x1+dn, x0+dn), s, 2*r, s, 2*r));
		append(new Arc2D(r,th,th+M_PI));
		append(new ParamStretcher(new Line2D(x0-dn, x1-dn), s, 2*r, s, 2*r));
		append(new Arc2D(r,th+M_PI,th+2*M_PI));
	} else {
		append(new ParamStretcher(new Line2D(x1-dn, x0-dn), s, 2*r, s, 2*r));
		append(new Arc2D(-r,th+M_PI,th));
		append(new ParamStretcher(new Line2D(x0+dn, x1+dn), s, 2*r, s, 2*r));
		append(new Arc2D(-r,th,th-M_PI));
	}
}



//
//
//

/*
/// Parameter distorter for speeding up or slowing down passage through segment in PathJoiner
class PathStretcher: public DVFunc1<2,double> {
public:
	/// constructor
	PathStretcher(DVFunc1<2,double>* f0, double a, bool dodelete = true): u(a), f(f0), delete_f(dodelete) { assert(f); }
	/// destructor
	virtual ~PathStretcher() { if(delete_f) delete f; }
	/// evaluate function
	virtual vec2 operator()(double x) const { return (*f)(distort(x)); }
	/// derivative
	virtual vec2 deriv(double x) const { return f->deriv(distort(x)) * d_distort(x); }

	double u;				//< distortion parameter, 1 for no distortion, 0 to infty for slow/fast stretch
	
protected:
	/// distortion map function
	virtual double distort(double l) const { return l * (l*(2*l-3)*(u-1) + u) ; }
	/// distortion derivative function
	virtual double d_distort(double l) const { return 6*l*(l-1)*(u-1) + u; }

	DVFunc1<2,double>* f;	//< function being distorted
	bool delete_f;			//< whether to delete function f on destruction
};
*/
