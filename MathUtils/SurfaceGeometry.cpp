#include "SurfaceGeometry.hh"
#include <cmath>

vec3 SurfaceGeometry::snorm(const vec2& p) const {
	return cross(deriv(p,0),deriv(p,1));
}

//--------------------------------------

vec3 CylSurfaceGeometry::operator()(const vec2& p) const {
	assert(fr);
	
	mdouble z = (fz ? (*fz)(p[0]) : p[0]);
	mdouble r = (*fr)(z);
	mdouble phi = 2*M_PI*p[1];
	mdouble c = cos(phi);
	mdouble s = sin(phi);
	
	return vec3(r*c,r*s,z);
}

vec3 CylSurfaceGeometry::deriv(const vec2& p, unsigned int i) const {
	assert(fr);
	
	mdouble z = (fz ? (*fz)(p[0]) : p[0]);
	mdouble phi = 2*M_PI*p[1];
	mdouble c = cos(phi);
	mdouble s = sin(phi);
	
	if(i==0) {
		mdouble dzdl = (fz ? fz->deriv(p[0]) : 1.);
		mdouble drdl = fr->deriv(z)*dzdl;
		return vec3(drdl*c,drdl*s,dzdl);
	}
	if(i==1) {
		mdouble r = (*fr)(z);
		return vec3(-r*s*2*M_PI, r*c*2*M_PI, 0);
	}
	
	assert(false);
	return vec3(0,0,0);
}

