#include "SurfaceProfiles.hh"
#include <cmath>

vec2 Arc2D::operator()(mdouble x) const {
	double th = t0+dt*x;
	return v + vec2(r*cos(th), r*sin(th));
}
vec2 Arc2D::deriv(mdouble x) const {
	double th = t0+dt*x;
	return vec2(-r*sin(th)*dt, r*cos(th)*dt);
}

//----

vec2 Dish::operator()(mdouble x) const {
	if(r>0) return vec2(z0-dz*((1-x)*(1-x)-1), r*(1-x));
	return vec2(z0-dz*(x*x-1), r*x);
}

//----

ParamStretcher::ParamStretcher(DVFunc1<2,mdouble>* f0, mdouble d0, mdouble r0, mdouble d1, mdouble r1, bool dodelete):
ri0(1./fabs(r0)), ri1(1./fabs(r1)), c0(-(d0-1)*fabs(r0)), c1((d1-1)*fabs(r1)), k(1), h(0), f(f0), delete_f(dodelete) {
	k = 1./(distort(1)-distort(0));
	h = -distort(0);
}

//----

RoundedSlab::RoundedSlab(mdouble z0, mdouble r, mdouble w): PathJoiner<2,mdouble>() {
	append(new ParamStretcher(new Line2D(vec2(z0-0.5*w, 0), vec2(z0-0.5*w, r-0.5*w)),
							  1, 1, 0.3, r));
	append(new Arc2D(0.5*w));
	append(new ParamStretcher(new Line2D(vec2(z0+0.5*w, r-0.5*w), vec2(z0+0.5*w, 0)),
							  0.3, r, 1, 1));
}

//----

RoundedTube::RoundedTube(vec2 x0, vec2 x1, mdouble r, mdouble endfrac): PathJoiner<2,mdouble>() {
	vec2 dl = (x1-x0).normalized()*fabs(r);
	vec2 dn = rhOrtho(dl);
	double th = atan2(dn[1],dn[0]);
	
	x1 -= dl;
	x0 += dl;
	mdouble sidelen = (x1-x0).mag();
	mdouble endlen = fabs(r)*M_PI;
	mdouble s = (1-endfrac)/endfrac * endlen/sidelen;
	
	if(r>0) {
		append(new ParamStretcher(new Line2D(x1+dn, x0+dn), s, r, s, r));
		append(new Arc2D(r,th,th+M_PI));
		append(new ParamStretcher(new Line2D(x0-dn, x1-dn), s, r, s, r));
		append(new Arc2D(r,th+M_PI,th+2*M_PI));
	} else {
		append(new ParamStretcher(new Line2D(x1-dn, x0-dn), s, r, s, r));
		append(new Arc2D(-r,th+M_PI,th));
		append(new ParamStretcher(new Line2D(x0+dn, x1+dn), s, r, s, r));
		append(new Arc2D(-r,th,th-M_PI));
	}
}



//
//
//

/*
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
*/
