#include "SurfaceGeometry.hh"
#include <cmath>

vec3 SurfaceGeometry::snorm(const vec2& p, bool normalized) const {
	if(normalized) {
		vec3 dx = deriv(p,0).normalized();
		vec3 dy = deriv(p,1).normalized();
		return cross(dx,dy);
	} else {
		return cross(deriv(p,0), deriv(p,0));
	}
}

Matrix<3,3,mdouble> SurfaceGeometry::rotToLocal(const vec2& x) const {
	vec3 v0 = deriv(x,0).normalized();
	vec3 v1 = deriv(x,1).normalized();
	vec3 v2 = cross(v0,v1);
	Matrix<3,3,mdouble> M;
	for(unsigned int i=0; i<3; i++) {
		M(0,i) = v0[i];
		M(1,i) = v1[i];
		M(2,i) = v2[i];
	}
	return M;
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

