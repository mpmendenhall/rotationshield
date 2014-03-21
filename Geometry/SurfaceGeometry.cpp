#include "SurfaceGeometry.hh"
#include "Integrator.hh"
#include <cmath>
#include <iostream>
#include <algorithm>

vec3 SurfaceGeometry::snorm(const vec2& p, bool normalized) const {
	if(normalized) {
		vec3 dx = deriv(p,0).normalized();
		vec3 dy = deriv(p,1).normalized();
		return cross(dx,dy);
	} else {
		return cross(deriv(p,0), deriv(p,1));
	}
}

vec2 SurfaceGeometry::d_pathlength(vec2 l) const {
	return vec2(deriv(l,0).mag(), deriv(l,1).mag());
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

mdouble SurfaceGeometry::dA(const vec2& l) const {
	return snorm(l,false).mag();
}

mdouble integration_dA(vec2 l, void* params) {
	SurfaceGeometry* S = (SurfaceGeometry*)params;
	return S->dA(l);
}

mdouble SurfaceGeometry::area(const vec2& ll, const vec2& ur) {
	return myIntegrator.integrate2D(&integration_dA, ll, ur, this);
}

void SurfaceGeometry::cache_sincos(double theta, double& s, double& c) const {
	static double c_th = theta;
	static double c_s = sin(theta);
	static double c_c = cos(theta);
	if(c_th != theta) { c_th=theta; c_c=cos(c_th); c_s=sin(c_th); }
	s = c_s;
	c = c_c;
}

void SurfaceGeometry::proximity(vec3 p, vec2 ll, vec2 ur, double& mn, double& mx) const {
	const int nx = 4;
	const int ny = 4;
	mdouble r[nx*ny];
	
	unsigned int i=0;
	for(int x = 0; x<nx; x++) {
		for(int y = 0; y<ny; y++) {
			double l1 = double(x)/(nx-1);
			double l2 = double(y)/(ny-1);
			vec2 v( ll[0]*(1-l1)+ur[0]*l1, ll[1]*(1-l2)+ur[1]*l2);
			r[i++] = ( p - (*this)(v) ).mag2();
		}
	}
	
	mn = sqrt(*std::max_element(r,r+nx*ny));
	mx = sqrt(*std::min_element(r,r+nx*ny));
}

//--------------------------------------

vec2 CylSurfaceGeometry::cache_profile(mdouble l) const {
	assert(zr_profile);
	static mdouble c_l = l;
	static vec2 c_v = (*zr_profile)(l);
	if(c_l != l) { c_l = l; c_v = (*zr_profile)(l); }
	return c_v;
}

vec3 CylSurfaceGeometry::operator()(const vec2& p) const {
	vec2 zr = cache_profile(p[0]);
	mdouble phi = 2*M_PI*p[1];
	double s,c;
	cache_sincos(phi,s,c);
	return vec3(zr[1]*c, zr[1]*s, zr[0]);
}

vec3 CylSurfaceGeometry::deriv(const vec2& p, unsigned int i) const {

	assert(zr_profile);
	vec2 zr = cache_profile(p[0]);
	mdouble phi = 2*M_PI*p[1];
	double s,c;
	cache_sincos(phi,s,c);
	
	if(i==0) {
		vec2 dzr = zr_profile->deriv(p[0]);
		assert(dzr.mag2());
		return vec3(dzr[1]*c, dzr[1]*s, dzr[0]);
	} else if(i==1) {
		return vec3(-zr[1]*s*2*M_PI, zr[1]*c*2*M_PI, 0);
	}
	
	assert(false);
	return vec3(0,0,0);
}

mdouble CylSurfaceGeometry::dA(const vec2& l) const {
	vec2 zr = cache_profile(l[0]);
	vec2 dzr = zr_profile->deriv(l[0]);
	return 2*M_PI*zr[1]*dzr.mag();
}

double integration_cyl_dA(double x, void* params) {
	CylSurfaceGeometry* S = (CylSurfaceGeometry*)params;
	return S->dA(vec2(x,0));
}

mdouble CylSurfaceGeometry::area(const vec2& ll, const vec2& ur) {
	return myIntegrator.integrate(&integration_cyl_dA, ll[0], ur[0], this)*(ur[1]-ll[1]);
}

//--------------------------------------------


