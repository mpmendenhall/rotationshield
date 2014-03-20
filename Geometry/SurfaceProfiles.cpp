#include "SurfaceProfiles.hh"

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

RoundedSlab::RoundedSlab(mdouble z0, mdouble r, mdouble w): PathJoiner<2,mdouble>() {
	append(new Line2D(vec2(z0-0.5*w, 0), vec2(z0-0.5*w, r-0.5*w)));
	append(new PathStretcher(new Arc2D(0.5*w), 2.5));
	append(new Line2D(vec2(z0+0.5*w, r-0.5*w), vec2(z0+0.5*w, 0)));
}

//----

