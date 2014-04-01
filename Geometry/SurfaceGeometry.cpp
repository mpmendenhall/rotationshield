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

Matrix<3,3,double> SurfaceGeometry::rotToLocal(const vec2& x) const {
	vec3 v0 = deriv(x,0).normalized();
	vec3 v1 = deriv(x,1).normalized();
	vec3 v2 = cross(v0,v1);
	Matrix<3,3,double> M;
	for(unsigned int i=0; i<3; i++) {
		M(0,i) = v0[i];
		M(1,i) = v1[i];
		M(2,i) = v2[i];
	}
	return M;
}

double SurfaceGeometry::dA(const vec2& l) const {
	return snorm(l,false).mag();
}

double integration_dA(mvec l, void* params) {
	SurfaceGeometry* S = (SurfaceGeometry*)params;
	return S->dA(vec2(l[0],l[1]));
}

double SurfaceGeometry::area(const vec2& ll, const vec2& ur) const {
	return myIntegratorND.integrate(&integration_dA, ll, ur, (void*)this);
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
	double r[nx*ny];
	
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

struct f2_to_fN_params {
	mvec (*f2)(vec2, void*);
	void* fparams;
};

mvec f2_to_fN(mvec v, void* pp) {
	f2_to_fN_params* p = (f2_to_fN_params*)pp;
	return (*p->f2)(vec2(v[0],v[1]), p->fparams);
}

mvec SurfaceGeometry::subdividedIntegral(mvec (*f)(vec2, void*), unsigned int fdim, void* fparams, vec2 ll, vec2 ur, unsigned int ndx, unsigned int ndy) const {
	
	ndx = ndx?ndx:dflt_integrator_ndivs_x;
	ndy = ndy?ndx:dflt_integrator_ndivs_y;
	for(unsigned int i=0; i<2; i++)
		ur[i] = (ur[i]-ll[i])/(i?ndy:ndx);

	mvec m;
	for(unsigned int nx=0; nx<ndx; nx++) {
		for(unsigned int ny=0; ny<ndy; ny++) {
			vec2 lll = ll + ur*vec2(nx,ny);
			mvec mi;
			if(polar_integral_center) {
				mi = myIntegratorND.integratePolar(f, fdim, *polar_integral_center, lll, lll+ur, fparams, -666, polar_r0);
			} else if(myIntegrator2D.getMethod() == INTEG_GSL_CQUAD) {
				mi = myIntegrator2D.integrate2D(f, lll, lll+ur, fparams);
			} else {
				f2_to_fN_params p;
				p.f2 = f;
				p.fparams = fparams;
				mi = myIntegratorND.integrate(&f2_to_fN, fdim, mvec(lll), mvec(lll+ur), &p);
			}
 			if(!nx && !ny) m = mi;
			else m += mi;
		}
	}
	return m;
}

//--------------------------------------

vec2 CylSurfaceGeometry::cache_profile(double l) const {
	assert(zr_profile);
	static double c_l = l;
	static vec2 c_v = (*zr_profile)(l);
	if(c_l != l) { c_l = l; c_v = (*zr_profile)(l); }
	return c_v;
}

vec3 CylSurfaceGeometry::operator()(const vec2& p) const {
	vec2 zr = cache_profile(p[0]);
	double phi = 2*M_PI*p[1];
	double s,c;
	cache_sincos(phi,s,c);
	return vec3(zr[1]*c, zr[1]*s, zr[0]);
}

vec3 CylSurfaceGeometry::deriv(const vec2& p, unsigned int i) const {

	assert(zr_profile);
	vec2 zr = cache_profile(p[0]);
	double phi = 2*M_PI*p[1];
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

double CylSurfaceGeometry::dA(const vec2& l) const {
	vec2 zr = cache_profile(l[0]);
	vec2 dzr = zr_profile->deriv(l[0]);
	return 2*M_PI*zr[1]*dzr.mag();
}

double integration_cyl_dA(double x, void* params) {
	CylSurfaceGeometry* S = (CylSurfaceGeometry*)params;
	return S->dA(vec2(x,0));
}

double CylSurfaceGeometry::area(const vec2& ll, const vec2& ur) const {
	return myIntegrator.integrate(&integration_cyl_dA, ll[0], ur[0], (void*)this)*(ur[1]-ll[1]);
}

//--------------------------------------------


