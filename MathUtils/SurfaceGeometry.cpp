#include "SurfaceGeometry.hh"
#include <cmath>

vec3 SurfaceGeometry::snorm(const vec2& p) const {
	return cross(deriv(p,0),deriv(p,1));
}

//--------------------------------------

vec3 CylSurfaceGeometry::operator()(const vec2& p) const {
	assert(fprofile);
	vec2 zr = (*fprofile)(p[0]);
	mdouble phi = 2*M_PI*p[1];
	mdouble c = cos(phi);
	mdouble s = sin(phi);
	return vec3(zr[1]*c, zr[1]*s, zr[0]);
}

vec3 CylSurfaceGeometry::deriv(const vec2& p, unsigned int i) const {

	assert(fprofile);
	vec2 zr = (*fprofile)(p[0]);
	mdouble phi = 2*M_PI*p[1];
	mdouble c = cos(phi);
	mdouble s = sin(phi);
	
	if(i==0) {
		vec2 dzr = fprofile->deriv(p[0]);
		return vec3(dzr[1]*c, dzr[1]*s, dzr[0]);
	}
	if(i==1) {
		return vec3(-zr[1]*s*2*M_PI, zr[1]*c*2*M_PI, 0);
	}
	
	assert(false);
	return vec3(0,0,0);
}

